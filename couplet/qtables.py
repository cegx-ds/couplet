#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
from collections import defaultdict
import logging
import math

import numpy as np
import pandas as pd
import pysam
import scipy.stats as stats
import tqdm

from couplet.rules import five_bp_rule, six_bp_rule

_five_bp_modC_codes = (("T", "C"), ("C", "C"))
_six_bp_modC_codes = (
    ("T", "C"),
    ("C", "C"),
    ("G", "A"),
    ("G", "G"),
)


def _reverse_rule(rule):
    return {value: key for key, value in rule.items() if value[0] != "N"}


def load_regions(regions_file, contig=None):
    regions_df = pd.read_csv(regions_file, sep="\t", header=None)

    # Bed files can have more than 3 columns, so subset and rename first 3 columns
    regions_df = regions_df[[0, 1, 2]]
    regions_df = regions_df.rename({0: "contig", 1: "start", 2: "end"}, axis="columns")

    if contig is not None:
        # Subset the .bed file by chromosome/contig of interest
        regions_df = regions_df[regions_df["contig"] == contig]
    return [(row.contig, row.start, row.end) for row in regions_df.itertuples()]


def _error_prob_to_phred(error_prob):

    # Temporarily suspend division by zero warnings. These can occur in np.log10(), but
    # are handled in np.nan_to_num(..., posinf=60).
    np.seterr(divide="ignore")
    phreds = np.rint(np.nan_to_num(-10 * np.log10(error_prob), posinf=60)).astype(int)
    np.seterr(divide="warn")

    return phreds


def make_qtable(
    bam_path,
    fasta_path,
    num_reads=1000000,
    rule="5bp",
    regions=None,
    contig=None,
    mapq_threshold=20,
    aggregate_mod_counts=True,
):
    """Estimate the qtable from empirical data

    Arguments:
        bam_path (str or pathlib.Path): path to the input BAM file
        fasta_path (str or pathlib.Path): path to the reference FASTA file
        rule (str): the name of the resolution rule
        regions (list): a list of (contig, start, end) tuples describing the regions
            to use to estimate the q-table
        contig (str): the name of the contig on which to estimate the q-table
        mapq_threshold (float): minimum mapping quality for a read to be included
        num_reads (int): how many reads to use to estimate the q-table
    """
    bam = pysam.AlignmentFile(bam_path, "rb", reference_filename=fasta_path)
    fasta = pysam.FastaFile(fasta_path)

    if rule == "5bp":
        rule_dict = _reverse_rule(five_bp_rule)

    elif rule == "6bp":
        rule_dict = _reverse_rule(six_bp_rule)

    else:
        assert False, f"{rule} is not a valid"

    correct = defaultdict(int)
    incorrect = defaultdict(int)

    if regions is None:
        if contig is not None:
            ref_id = bam.references.index(contig)
            regions = [(contig, 0, bam.lengths[ref_id])]
        else:
            regions = [
                (contig, 0, end) for contig, end in zip(bam.references, bam.lengths)
            ]

    processed = 0
    for contig, start, end in tqdm.tqdm(regions):
        ref = fasta.fetch(contig, start, end)
        for read in bam.fetch(contig=contig, start=start, stop=end):
            if read.mapping_quality < mapq_threshold:
                continue

            query_sequence = read.query_sequence

            if read.is_reverse:
                # Get the query sequence in the original fastq orientation
                query_sequence_rev_comp = read.get_forward_sequence()

            first_q = pysam.qualitystring_to_array(read.get_tag("XQ"))
            second_q = pysam.qualitystring_to_array(read.get_tag("YQ"))
            epi_seq = read.get_tag("XE")

            for query_pos, ref_pos in read.get_aligned_pairs():

                if (
                    query_pos is None
                    or ref_pos is None
                    or ref_pos < start
                    or ref_pos >= end
                    or ref[ref_pos - start] == "N"
                    or query_sequence[query_pos] == "N"
                ):
                    continue

                if read.is_reverse:
                    first_base, second_base = rule_dict[
                        (
                            query_sequence_rev_comp[-(query_pos + 1)],
                            epi_seq[-(query_pos + 1)],
                        )
                    ]
                    key = (
                        first_base,
                        second_base,
                        first_q[-(query_pos + 1)],
                        second_q[-(query_pos + 1)],
                    )

                else:
                    first_base, second_base = rule_dict[
                        (query_sequence[query_pos], epi_seq[query_pos])
                    ]
                    key = (
                        first_base,
                        second_base,
                        first_q[query_pos],
                        second_q[query_pos],
                    )

                if query_sequence[query_pos] == ref[ref_pos - start]:
                    correct[key] += 1
                else:
                    incorrect[key] += 1
            processed += 1
            if processed >= num_reads:
                break
        logging.debug("Processed %d reads.", processed)
        if processed >= num_reads:
            logging.info("Processed %d reads. Stopping", processed)
            break
    logging.info("Estimated Q-table from %d reads", processed)

    # If requested, combine the correct/incorrect counts for C & modC/ G * mC and compute
    # common counts.
    if aggregate_mod_counts:
        correct, incorrect = _aggregate_mod_counts(correct, incorrect, rule)

    qtable = pd.DataFrame(
        data={
            "first_base": pd.Series(dtype=str),
            "second_base": pd.Series(dtype=str),
            "first_phred": pd.Series(dtype=int),
            "second_phred": pd.Series(dtype=int),
            "correct": pd.Series(dtype=int),
            "incorrect": pd.Series(dtype=int),
        }
    )

    for i, key in enumerate(set(correct) | set(incorrect)):
        first_base, second_base, first_phred, second_phred = key
        num_correct = correct[key]
        num_incorrect = incorrect[key]
        qtable.loc[i] = [
            first_base,
            second_base,
            first_phred,
            second_phred,
            num_correct,
            num_incorrect,
        ]
    qtable["error_rate"] = qtable.incorrect / (qtable.correct + qtable.incorrect)
    qtable["phred"] = _error_prob_to_phred(qtable["error_rate"])

    return qtable, processed


def aggregate_qtables(qtable_list):
    """Build an aggregate quality table from a list of input quality tables"""
    reindexed_tables = [
        df.set_index(["first_base", "second_base", "first_phred", "second_phred"])
        for df in qtable_list
    ]
    combined_table = reindexed_tables[0][["correct", "incorrect"]]
    for df in reindexed_tables[1:]:
        combined_table = combined_table.add(df[["correct", "incorrect"]], fill_value=0)
        print(combined_table.correct)
    result = combined_table.reset_index()
    result["error_rate"] = result.incorrect / (result.correct + result.incorrect)
    result["phred"] = _error_prob_to_phred(result["error_rate"])

    return result


def _aggregate_mod_counts(correct, incorrect, rule):
    """
    This function takes correct and incorrect count dictionaries and aggregate the
    counts for C and modC bases in 5L mode, and also for both Gs (G or mCG) in 6L mode.
    That is, for each of correct/incorrect dictionaries, the
      - (C and modC) or (G and mC) counts are combined depending on 5L/6L mode
      - Other counts are left unchanged
      - Counts unique to (i.e.) C/modC are made "shared", e.g. if C has a count of 1292 before
        aggregation and modC has no counts before aggregation, then both C and modC have
        counts of 1292 after aggregation.
      - Counts unique to correct/incorrect dictionaries stay unique to this dictionary

    Args:
      - correct: defaultdict of correct counts. Of form: { ["T","C",37,24]: 1292, ...}
      - incorrect: defaultdict of incorrect counts.
      - rule: (str) the resolution rule.

    Returns:
       correct and incorrect dictionaries with counts aggregated across C and mod C.
    """

    if rule == "5bp":
        modC_codes = _five_bp_modC_codes

    elif rule == "6bp":
        modC_codes = _six_bp_modC_codes

    else:
        assert False, f"{rule} is not a valid"

    # Get dictionaries with only C/modC entries
    correct_modC = {k: v for k, v in correct.items() if k[0:2] in modC_codes}
    incorrect_modC = {k: v for k, v in incorrect.items() if k[0:2] in modC_codes}

    # Create a "shared" dictionary with summed counts. Counts are summed for all entries
    # with identical (Q1, Q2).
    G_mC_shared_correct = defaultdict(int)
    C_modC_shared_correct = defaultdict(int)

    # Take correct R1/R2 pairings
    for k, v in correct_modC.items():

        # If pairing specific to 6L (i.G/mC)
        if k[0:2] in _six_bp_modC_codes and k[0:2] not in _five_bp_modC_codes:

            # Create shared dict for mC pairings
            G_mC_shared_correct[("G_mC", "G_mC") + k[2:4]] += v

        # Else pairing is for C/modC (5L)
        else:
            C_modC_shared_correct[("C_modC", "C_modC") + k[2:4]] += v

    # Init dicts
    G_mC_shared_incorrect = defaultdict(int)
    C_modC_shared_incorrect = defaultdict(int)

    # Take incorrect R1/R2 pairings
    for k, v in incorrect_modC.items():

        # If pairing specific to 6L (i.G/mC)
        if k[0:2] in _six_bp_modC_codes and k[0:2] not in _five_bp_modC_codes:

            # Group Phreds together and assign both to G_mC
            G_mC_shared_incorrect[("G_mC", "G_mC") + k[2:4]] += v

        # Else pairing is for C/modC (5L)
        else:
            C_modC_shared_incorrect[("C_modC", "C_modC") + k[2:4]] += v

    # Update original dictionaries: Each pair of C/modC codes get the same value for a
    # given (Q1, Q2) entry.
    # Loop over modC pairings
    for pair in modC_codes:

        # Check if 6l/(G/mC)
        if pair in _six_bp_modC_codes and pair not in _five_bp_modC_codes:

            # Update counts with combined counts
            for k, v in G_mC_shared_correct.items():
                correct[pair + k[2:4]] = v
            for k, v in G_mC_shared_incorrect.items():
                incorrect[pair + k[2:4]] = v

        # Repeat for 5L modC
        else:
            for k, v in C_modC_shared_correct.items():
                correct[pair + k[2:4]] = v

            for k, v in C_modC_shared_incorrect.items():
                incorrect[pair + k[2:4]] = v

    return correct, incorrect
