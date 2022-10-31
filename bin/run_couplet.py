#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
"""
The script resolves two SingleShot reads
Usage: python run_couplet.py  --fq1=TEST_R1.fq.gz  --fq2=TEST_R2.fq.gz --mismatch_threshold=0.05 --phred=prob

"""

import argparse
import logging
import re
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import sys
import os
import functools

from couplet.utils import (
    get_base_name,
    get_outfile_names,
    check_python_version,
    load_quality_table,
)
from couplet.export import (
    log_core_stats,
    log_additional_stats,
    log_and_plot_additional_stats_single,
)
from couplet.resolve import (
    resolve_phred_min,
    resolve_phred_prob,
    resolve_phred_with_qtable,
)
from couplet.rules import (
    five_bp_rule,
    five_bp_error_codes,
    five_bp_modifications,
    six_bp_rule,
    six_bp_error_codes,
    six_bp_modifications,
    ResolutionRule,
)
from couplet.core import resolve_read_pair, update_stats


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fq1",
        help="forward read names",
    )

    parser.add_argument(
        "--fq2",
        help="reverse read names",
    )

    parser.add_argument(
        "--mismatch-threshold",
        help="Enter a value between 0 and 1 for allowed fraction of mismatches",
        default=0.05,
        type=float,
    )
    parser.add_argument(
        "--phred",
        help="Enter a mode: min or prob",
        default="prob",
        choices=["min", "prob", "qtable"],
        type=str,
    )
    parser.add_argument(
        "--quality-table",
        help="Enter the path to a quality table",
        default=None,
        type=str,
    )

    parser.add_argument(
        "--rule",
        help="Enter a rule name: 5bp or 6bp",
        default="5bp",
        choices=["5bp", "6bp"],
        type=str,
    )

    parser.add_argument(
        "--gap-penalty",
        default=-2,
        help="Gap penalty parameter for NW alignment default is -2",
        type=float,
    )
    parser.add_argument(
        "--end-gap-penalty",
        default=-5,
        help="gap penalty for begining and end gaps  for NW alignment default is -5",
        type=float,
    )
    parser.add_argument(
        "--match-award",
        default=1,
        help="Match award parameter for NW alignment default is 1",
        type=float,
    )
    parser.add_argument(
        "--mismatch-penalty",
        default=-1,
        help="Mismatch  penalty parameter for NW alignment default is -1",
        type=float,
    )
    parser.add_argument(
        "--no-mismatch-aware-trimming",
        dest="mismatch_aware_trimming",
        action="store_false",
        help="Flag to turn off mismatch aware trimming.",
    )
    parser.add_argument(
        "--min-read-length",
        default=50,
        help="Minimum read length required to not discard a read.",
        type=float,
    )
    parser.add_argument(
        "--additional-stats",
        dest="additional_stats",
        action="store_true",
        help="Flag to turn on logging of additional statistics. This the default behaviour.",
    )
    parser.add_argument(
        "--no-additional-stats",
        dest="additional_stats",
        action="store_false",
        help="Flag to turn off logging of additional statistics.",
    )
    parser.add_argument(
        "--generate-plots",
        dest="generate_plots",
        action="store_true",
        help="Flag to indicate that plotting should take place here, rather than in post-processing.",
    )
    parser.add_argument(
        "--episeq-in-qname",
        dest="episeq_in_qname",
        action="store_true",
        help="Flag to export the epigenetic sequence as part of the fastq read name (or QNAME in .sam specification).",
    )
    parser.add_argument(
        "--orig-quals-in-sam-tag",
        dest="orig_quals_in_sam_tag",
        action="store_true",
        help="Flag to turn on storing of the original quality scores (Q1, Q2) in a SAM tag (XQ and YQ respectively).",
    )

    parser.add_argument(
        "--xe-tag",
        dest="xe_tag",
        action="store_true",
        help="Flag to turn on storing of modification information in a CEGX tag (XE).",
    )

    parser.set_defaults(
        episeq_in_qname=False,
        generate_plots=False,
        additional_stats=True,
        mismatch_aware_trimming=True,
        xe_tag=False,
        orig_quals_in_sam_tag=False,
    )

    args = parser.parse_args()
    if len(sys.argv) == 1:
        # display help message when no args are passed.
        parser.print_help()
        sys.exit()

    # Check Python version is above 3.9.0
    check_python_version([3, 9, 0])

    # Setup logging
    log_output_file = get_base_name(args) + "_couplet.log"
    logging.basicConfig(
        filename=log_output_file,
        format="%(asctime)-15s %(filename)s:%(lineno)d %(message)s",
    )
    LOGGER = logging.getLogger("root")
    LOGGER.setLevel(logging.INFO)
    LOGGER.addHandler(logging.StreamHandler())

    # Log all input arguments
    LOGGER.info("Command line inputs:")
    for k, v in vars(args).items():
        LOGGER.info(k + ": " + str(v))

    # Choosing mode phred calc
    if args.phred == "min":
        resolve_phred_fn = resolve_phred_min
    elif args.phred == "qtable":
        if args.quality_table is None:
            LOGGER.error(
                "Requested 'qtable' mode for Phred score resolution, but no quality table was provided. Exiting."
            )
            sys.exit(1)
        elif os.path.exists(args.quality_table):
            LOGGER.info("Using quality table: " + args.quality_table)
            quality_table = load_quality_table(args.quality_table)
            resolve_phred_fn = functools.partial(
                resolve_phred_with_qtable, quality_table=quality_table
            )
        else:
            LOGGER.error(
                f"Requested 'qtable' mode for Phred score resolution, but quality table '{args.quality_table}' does not exist. Exiting."
            )
            sys.exit(1)
    else:
        resolve_phred_fn = resolve_phred_prob

    if args.rule == "5bp":
        rule = ResolutionRule(
            five_bp_rule,
            five_bp_error_codes,
            five_bp_modifications,
            match_award=args.match_award,
            gap_penalty=args.gap_penalty,
            mismatch_penalty=args.mismatch_penalty,
            end_gap_penalty=args.end_gap_penalty,
        )
    elif args.rule == "6bp":
        rule = ResolutionRule(
            six_bp_rule,
            six_bp_error_codes,
            six_bp_modifications,
            match_award=args.match_award,
            gap_penalty=args.gap_penalty,
            mismatch_penalty=args.mismatch_penalty,
            end_gap_penalty=args.end_gap_penalty,
        )
    else:
        print("Not defined !")
        sys.exit()

    (
        out_discard_file1,
        out_discard_file2,
        resolved_path,
        stats_file,
        additional_stats_file,
    ) = get_outfile_names(args)

    LOGGER.info(f"Writing resolved read to: {resolved_path}")
    stats = {}
    with gzip.open(out_discard_file1, mode="wt") as out_discard_file1F, gzip.open(
        out_discard_file2, mode="wt"
    ) as out_discard_file2F, gzip.open(resolved_path, mode="wt") as resolved_handle:
        read1_iter = SeqIO.parse(gzip.open(args.fq1, "rt"), "fastq")
        read2_iter = SeqIO.parse(gzip.open(args.fq2, "rt"), "fastq")
        for r1, r2 in zip(read1_iter, read2_iter):
            result, read_stats = resolve_read_pair(
                r1,
                r2,
                rule,
                resolve_phred_fn=resolve_phred_fn,
                mismatch_threshold=args.mismatch_threshold,
                use_mismatch_aware_trimming=args.mismatch_aware_trimming,
                min_read_length=args.min_read_length,
                episeq_in_qname=args.episeq_in_qname,
                log_additional_stats=args.additional_stats,
                orig_quals_in_sam_tag=args.orig_quals_in_sam_tag,
                xe_tag=args.xe_tag,
            )
            read_stats["num_reads"] = 1
            stats = update_stats(stats, read_stats)
            if result is not None:
                SeqIO.write(result, resolved_handle, "fastq")
            else:
                SeqIO.write(r1, out_discard_file1F, "fastq")
                SeqIO.write(r2, out_discard_file2F, "fastq")

    # Log and plot statistics
    log_core_stats(stats, stats_file)
    if args.additional_stats:
        if args.generate_plots:
            log_and_plot_additional_stats_single(stats, rule, additional_stats_file)
        else:
            log_additional_stats(stats, rule, additional_stats_file)


if __name__ == "__main__":
    main()
