#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import argparse
import logging
import os

import pandas as pd

from couplet.qtables import load_regions, make_qtable


def main():
    parser = argparse.ArgumentParser("Estimate a Q-table for phred score qtables")
    parser.add_argument("--bam-path", help="path to the BAM file", required=True)
    parser.add_argument(
        "--ref-path", help="path to the reference FASTA file", required=True
    )
    parser.add_argument(
        "--pipeline-version", help="Pipeline version (for provenance)", required=True
    )
    parser.add_argument(
        "--output-dir", help="directory in which to write any output files", default="."
    )
    parser.add_argument(
        "--regions-file", help="path to a file containing selected regions"
    )
    parser.add_argument(
        "--contig", help="list of contigs from which to estimate the Q-table"
    )
    parser.add_argument("--rule", help="which SingleShot rule to use", default="5bp")
    parser.add_argument(
        "--tag", help="a tag to include in the output file name", default=""
    )
    parser.add_argument(
        "--mapping-quality",
        help="minimum mapping quality for a read to be included",
        default=30,
        type=int,
    )
    parser.add_argument(
        "--max-reads",
        help="maximum number of reads to process",
        default=1000000,
        type=int,
    )
    parser.add_argument(
        "--separate-mod-counts",
        dest="aggregate_mod_counts",
        action="store_false",
        help="flag to turn off aggregation of counts for modifications of the same base",
    )
    parser.set_defaults(
        aggregate_mod_counts=True,
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    if args.contig is not None:
        tag = f"{args.tag}_{args.contig}_{args.rule}_q{args.mapping_quality}"
    else:
        tag = f"{args.tag}_{args.rule}_q{args.mapping_quality}"
    log_file = os.path.join(args.output_dir, f"{tag}_make_qtable.log")
    logging.basicConfig(
        filename=log_file,
        format="%(asctime)-15s %(filename)s:%(lineno)d %(message)s",
        level=logging.INFO,
    )
    logging.getLogger().addHandler(logging.StreamHandler())

    if args.regions_file is not None:
        logging.info("Using high confidence regions: %s", args.regions_file)
        regions = load_regions(args.regions_file, args.contig)
    else:
        regions = None

    qtable, num_reads = make_qtable(
        args.bam_path,
        args.ref_path,
        rule=args.rule,
        regions=regions,
        contig=args.contig,
        mapq_threshold=args.mapping_quality,
        num_reads=args.max_reads,
        aggregate_mod_counts=args.aggregate_mod_counts,
    )

    qtable = qtable.sort_values(
        by=["first_phred", "second_phred", "first_base", "second_base"], ascending=False
    )

    output_file = os.path.join(args.output_dir, f"{tag}_qtable.csv")
    with open(output_file, "w") as out:
        out.write(f"# Pipeline version: {args.pipeline_version}\n")
        out.write(f"# Input file: {os.path.basename(args.bam_path)}\n")
        if args.regions_file is not None:
            out.write(f"# Regions file: {os.path.basename(args.regions_file)}\n")
        if args.contig is not None:
            out.write(f"# Contig: {args.contig}\n")
        out.write(f"# Rule: {args.rule}\n")
        out.write(f"# MAPQ threshold: {args.mapping_quality}\n")
        out.write(f"# Aggregate C / modC: {args.aggregate_mod_counts}\n")
        out.write(f"# Number of reads: {num_reads}\n")
    qtable.to_csv(output_file, index=False, mode="a")


if __name__ == "__main__":
    main()
