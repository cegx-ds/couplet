#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import csv
import pandas as pd
import pytest
import pysam
import numpy as np

from ssresolve.qtables import (
    aggregate_qtables,
    load_regions,
    make_qtable,
    _aggregate_mod_counts,
    _error_prob_to_phred,
)


def make_read(
    name,
    flag,
    ref_id,
    ref_start,
    mapq,
    cigar,
    distance,
    next_ref_id,
    next_ref_start,
    sequence,
    epi_seq,
    first_phred,
    second_phred,
):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    a.reference_id = ref_id
    a.reference_start = ref_start
    a.mapping_quality = mapq
    a.cigarstring = cigar
    a.next_reference_id = next_ref_id
    a.next_reference_start = next_ref_start
    # a.template_length = 167
    a.query_sequence = sequence
    a.query_qualities = pysam.qualitystring_to_array("<" * len(sequence))
    a.tags = (
        ("NM", distance),
        ("RG", "L1"),
        ("XE", epi_seq),
        ("XQ", first_phred),
        ("YQ", second_phred),
    )
    return a


# A small genome
@pytest.fixture
def fasta_path(tmp_path):
    fasta_path = tmp_path / "tmp.fasta"
    with open(fasta_path, "w") as out:
        out.write(">chr1\nAGCTNCG\n")
    pysam.faidx(fasta_path.as_posix())
    return fasta_path


# A set of reads covering the genome
@pytest.fixture
def bam_path(tmp_path):
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 7, "SN": "chr1"}]}

    bam_path = tmp_path / "tmp.bam"
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # Reads 1 to 4 contribute the bulk of the error rate stats
        # Genome: AGCTNPG
        # Read 1: AGCT
        read_1 = make_read(
            "read_1", 0, 0, 0, 20, "4M", 0, 0, 1, "AGCT", "AGCT", "AAAA", "AAAA"
        )
        outf.write(read_1)

        # Read 2 is a reverse read
        # Genome (reverse)             PGNAGCT
        # Read 2 (FASTQ orientation):    GAGT
        # Genome (forward)             AGCTNPG
        # Read 2 (BAM orientation):     ACTC
        read_2 = make_read(
            "read_2", 16, 0, 1, 20, "4M", 1, 0, 2, "ACTC", "GAGT", "AAAA", "AAAA"
        )
        outf.write(read_2)

        # Genome: AGCTNPG
        # Read 3:   PCGP
        read_3 = make_read(
            "read_3", 0, 0, 2, 20, "4M", 1, 0, 3, "CCGC", "PCGP", "AAAA", "AAAA"
        )
        outf.write(read_3)

        # Read 4 is a reverse read
        # Genome (reverse)             PGNAGCT
        # Read 2 (FASTQ orientation):  PGAC
        # Genome (forward)             AGCTNPG
        # Read 2 (BAM orientation):       GTCG
        # Genome assumes that CpG is methylated on both strands
        read_4 = make_read(
            "read_4", 16, 0, 3, 20, "4M", 1, 0, 4, "GTCG", "PGAC", "AAAA", "AAAA"
        )
        outf.write(read_4)

        # Read 5 and 6 check that indels and Ns are handled
        # Genome: AGCTNPG
        # Read 5:     T GPG
        read_5 = make_read(
            "read_5", 0, 0, 4, 20, "1M1D1M2S", 3, 0, 5, "TGCG", "TGPG", "AAAA", "AAAA"
        )
        outf.write(read_5)

        # Genome: AGCTNP-G
        # Read 6:      NCGG
        read_6 = make_read(
            "read_6", 0, 0, 5, 20, "1M1I1M1S", 3, 0, 6, "NCGG", "NCGG", "AAAA", "AAAA"
        )
        outf.write(read_6)

        # MAPQ < 20 (default), i.e. will be ignored
        # Genome: AGCTNPG
        # Read 7:       AAAA
        read_7 = make_read(
            "read_7", 0, 0, 6, 10, "1M3S", 3, 0, 7, "AAAA", "AAAA", "AAAA", "AAAA"
        )
        outf.write(read_7)

    pysam.index(bam_path.as_posix())
    return bam_path


# The reads in the BAM give the following pileups
# Bases from reads are arranged in columns
#
# Forward reads:
# 0 A A
# 1 G G
# 2 C CP
# 3 T TC
# 4 N  GT
# 5 P  P N
# 6 G   GG
#
# correct    A: 1, P: 2, G: 3, C: 1, T: 1
# incorrect: A: 0, P: 0, G: 0, C: 1, T: 0
#
# Reverse reads:
# 0 T
# 1 C T
# 2 G G
# 3 A A C
# 4 N G A
# 5 G   G
# 6 P   P
#
# correct:   A: 1, P: 1, G: 2, C: 0, T: 0
# incorrect: A: 0, P: 0, G: 0, C: 1, T: 1
#
# total correct:   A: 2, P: 3, G: 5, C: 1, T: 1
# total incorrect: A: 0, P: 0, G: 0, C: 2, T: 1
#
# These give the following empirical qtable
@pytest.fixture
def expected_qtable():
    qtable = pd.DataFrame(
        data=[
            ["A", "A", 32, 32, 2, 0, 0.0],  # Resolved A
            ["C", "C", 32, 32, 3, 0, 0.0],  # Resolved P
            ["G", "A", 32, 32, 5, 0, 0.0],  # Resolved G
            ["T", "C", 32, 32, 1, 2, 2 / 3],  # Resolved C
            ["T", "T", 32, 32, 1, 1, 1 / 2],  # Resolved T
        ],
        columns=[
            "first_base",
            "second_base",
            "first_phred",
            "second_phred",
            "correct",
            "incorrect",
            "error_rate",
        ],
    )
    return qtable


def test_makeq_table(bam_path, fasta_path, expected_qtable):
    sort_cols = ["first_base", "second_base", "first_phred", "second_phred"]
    actual_table, _ = make_qtable(bam_path, fasta_path, aggregate_mod_counts=False)
    actual_table = actual_table.sort_values(by=sort_cols).reset_index(drop=True)

    pd.testing.assert_frame_equal(
        actual_table[expected_qtable.columns],
        expected_qtable.sort_values(by=sort_cols).reset_index(drop=True),
    )


def test_makeq_table_set_contig(bam_path, fasta_path, expected_qtable):
    sort_cols = ["first_base", "second_base", "first_phred", "second_phred"]
    actual_table, _ = make_qtable(
        bam_path, fasta_path, contig="chr1", aggregate_mod_counts=False
    )
    actual_table = actual_table.sort_values(by=sort_cols).reset_index(drop=True)

    pd.testing.assert_frame_equal(
        actual_table[expected_qtable.columns],
        expected_qtable.sort_values(by=sort_cols).reset_index(drop=True),
    )


def test_makeq_table_set_region(bam_path, fasta_path, expected_qtable):
    sort_cols = ["first_base", "second_base", "first_phred", "second_phred"]
    actual_table, _ = make_qtable(
        bam_path, fasta_path, regions=[("chr1", 0, 7)], aggregate_mod_counts=False
    )
    actual_table = actual_table.sort_values(by=sort_cols).reset_index(drop=True)

    pd.testing.assert_frame_equal(
        actual_table[expected_qtable.columns],
        expected_qtable.sort_values(by=sort_cols).reset_index(drop=True),
    )


@pytest.fixture
def regions_path(tmp_path):
    regions_path = tmp_path / "regions.bed"
    with open(regions_path, "w") as fp:
        writer = csv.writer(fp, delimiter="\t")
        writer.writerow(("chr1", 0, 100))
        writer.writerow(("chr1", 500, 600))
        writer.writerow(("chr2", 0, 1000))
    return regions_path


def test_load_regions(regions_path):
    actual_regions = load_regions(regions_path)
    expected_regions = [("chr1", 0, 100), ("chr1", 500, 600), ("chr2", 0, 1000)]
    assert actual_regions == expected_regions


def test_load_regions_in_contig(regions_path):
    actual_regions = load_regions(regions_path, contig="chr1")
    expected_regions = [("chr1", 0, 100), ("chr1", 500, 600)]
    assert actual_regions == expected_regions


def test_aggregate_qtables(expected_qtable):
    input_qtables = [
        pd.DataFrame(
            data=[
                ["A", "A", 32, 32, 2, 0, 0.0, 60],  # Resolved A
                ["C", "C", 32, 32, 3, 0, 0.0, 60],  # Resolved P
                ["G", "A", 32, 32, 5, 0, 0.0, 60],  # Resolved G
                ["T", "C", 32, 32, 1, 2, 2 / 3, 1],  # Resolved C
            ],
            columns=[
                "first_base",
                "second_base",
                "first_phred",
                "second_phred",
                "correct",
                "incorrect",
                "error_rate",
                "phred",
            ],
        ),
        pd.DataFrame(
            data=[
                ["A", "A", 32, 32, 1, 1, 1 / 2, 3],  # Resolved A
                ["C", "C", 32, 32, 3, 1, 1 / 4, 6],  # Resolved P
                ["G", "A", 32, 32, 5, 3, 3 / 8, 4],  # Resolved G
                ["T", "C", 32, 32, 1, 1, 1 / 2, 3],  # Resolved C
                ["T", "T", 32, 32, 1, 1, 1 / 2, 3],  # Resolved T (only in second table)
            ],
            columns=[
                "first_base",
                "second_base",
                "first_phred",
                "second_phred",
                "correct",
                "incorrect",
                "error_rate",
                "phred",
            ],
        ),
    ]

    expected = pd.DataFrame(
        data=[
            ["A", "A", 32, 32, 3, 1, 1 / 4, 6],  # Resolved A
            ["C", "C", 32, 32, 6, 1, 1 / 7, 8],  # Resolved P
            ["G", "A", 32, 32, 10, 3, 3 / 13, 6],  # Resolved G
            ["T", "C", 32, 32, 2, 3, 3 / 5, 2],  # Resolved C
            ["T", "T", 32, 32, 1, 1, 1 / 2, 3],  # Resolved T
        ],
        columns=[
            "first_base",
            "second_base",
            "first_phred",
            "second_phred",
            "correct",
            "incorrect",
            "error_rate",
            "phred",
        ],
    )

    actual = aggregate_qtables(input_qtables)
    pd.testing.assert_frame_equal(actual, expected, check_dtype=False)


def test_aggregate_mod_counts_5bp():

    # This test tests that:
    #  - T/C and C/C counts are combined
    #  - Other counts are unchanged
    #  - Counts unique to T/C or C/C are made "shared"
    #  - Counts unique to correct/incorrect stays unique to the dictionary of relevance

    # Create input correct and incorrect count dictionaries
    correct = {}
    incorrect = {}

    correct[("A", "A", 37, 37)] = 8673
    correct[("A", "A", 24, 37)] = 4632
    correct[("C", "C", 37, 37)] = 1112
    correct[("C", "C", 24, 37)] = 877
    correct[("T", "C", 37, 37)] = 12315
    correct[("T", "C", 24, 37)] = 3026
    correct[("T", "C", 11, 24)] = 2881  # Unique to T/C and to correct

    incorrect[("A", "A", 37, 37)] = 8
    incorrect[("A", "A", 24, 37)] = 40
    incorrect[("C", "C", 37, 37)] = 13
    incorrect[("C", "C", 24, 37)] = 72
    incorrect[("C", "C", 11, 2)] = 7  # Unique to C/C and to incorrect
    incorrect[("T", "C", 37, 37)] = 12
    incorrect[("T", "C", 24, 37)] = 25

    # Define expectations
    correct_expectation = {}
    incorrect_expectation = {}

    correct_expectation[("A", "A", 37, 37)] = 8673
    correct_expectation[("A", "A", 24, 37)] = 4632
    correct_expectation[("C", "C", 37, 37)] = 13427
    correct_expectation[("C", "C", 24, 37)] = 3903
    correct_expectation[("C", "C", 11, 24)] = 2881
    correct_expectation[("T", "C", 37, 37)] = 13427
    correct_expectation[("T", "C", 24, 37)] = 3903
    correct_expectation[("T", "C", 11, 24)] = 2881

    incorrect_expectation[("A", "A", 37, 37)] = 8
    incorrect_expectation[("A", "A", 24, 37)] = 40
    incorrect_expectation[("C", "C", 37, 37)] = 25
    incorrect_expectation[("C", "C", 24, 37)] = 97
    incorrect_expectation[("C", "C", 11, 2)] = 7
    incorrect_expectation[("T", "C", 37, 37)] = 25
    incorrect_expectation[("T", "C", 24, 37)] = 97
    incorrect_expectation[("T", "C", 11, 2)] = 7

    # Run code and test
    correct_result, incorrect_result = _aggregate_mod_counts(
        correct, incorrect, rule="5bp"
    )

    assert correct_result == correct_expectation
    assert incorrect_result == incorrect_expectation


def test_aggregate_mod_counts_6bp():

    # This test tests that:
    #  - T/C and C/C counts are combined, and as are  G/G and G/A
    #  - Other counts are unchanged
    #  - Counts unique to above pairings (of pairs) are made "shared"
    #  - Counts unique to correct/incorrect stays unique to the dictionary of relevance

    # Create input correct and incorrect count dictionaries
    correct = {}
    incorrect = {}

    correct[("A", "A", 37, 37)] = 8673
    correct[("A", "A", 24, 37)] = 4632
    correct[("C", "C", 37, 37)] = 1112
    correct[("C", "C", 24, 37)] = 877
    correct[("T", "C", 37, 37)] = 12315
    correct[("T", "C", 24, 37)] = 3026
    correct[("T", "C", 11, 24)] = 2881  # Unique to T/C and to correct

    incorrect[("A", "A", 37, 37)] = 8
    incorrect[("A", "A", 24, 37)] = 40
    incorrect[("C", "C", 37, 37)] = 13
    incorrect[("C", "C", 24, 37)] = 72
    incorrect[("C", "C", 11, 2)] = 7  # Unique to C/C and to incorrect
    incorrect[("T", "C", 37, 37)] = 12
    incorrect[("T", "C", 24, 37)] = 25

    correct[("G", "A", 11, 24)] = 3000
    correct[("G", "A", 24, 37)] = 1789
    correct[("G", "A", 37, 37)] = 1350
    correct[("G", "G", 11, 24)] = 2006
    correct[("G", "G", 24, 37)] = 1995
    correct[("G", "G", 37, 37)] = 696
    correct[("G", "G", 37, 2)] = 800  # Unique to G/G and to correct

    incorrect[("G", "A", 11, 24)] = 10
    incorrect[("G", "A", 24, 37)] = 10
    incorrect[("G", "A", 37, 37)] = 10
    incorrect[("G", "A", 2, 11)] = 10  # Unique to G/A and to incorrect
    incorrect[("G", "G", 11, 24)] = 10
    incorrect[("G", "G", 24, 37)] = 10
    incorrect[("G", "G", 37, 37)] = 10

    # Define expectations
    correct_expectation = {}
    incorrect_expectation = {}

    correct_expectation[("A", "A", 37, 37)] = 8673
    correct_expectation[("A", "A", 24, 37)] = 4632
    correct_expectation[("C", "C", 37, 37)] = 13427
    correct_expectation[("C", "C", 24, 37)] = 3903
    correct_expectation[("C", "C", 11, 24)] = 2881
    correct_expectation[("T", "C", 37, 37)] = 13427
    correct_expectation[("T", "C", 24, 37)] = 3903
    correct_expectation[("T", "C", 11, 24)] = 2881

    incorrect_expectation[("A", "A", 37, 37)] = 8
    incorrect_expectation[("A", "A", 24, 37)] = 40
    incorrect_expectation[("C", "C", 37, 37)] = 25
    incorrect_expectation[("C", "C", 24, 37)] = 97
    incorrect_expectation[("C", "C", 11, 2)] = 7
    incorrect_expectation[("T", "C", 37, 37)] = 25
    incorrect_expectation[("T", "C", 24, 37)] = 97
    incorrect_expectation[("T", "C", 11, 2)] = 7

    correct_expectation[("G", "A", 11, 24)] = 5006
    correct_expectation[("G", "A", 24, 37)] = 3784
    correct_expectation[("G", "A", 37, 37)] = 2046
    correct_expectation[("G", "A", 37, 2)] = 800
    correct_expectation[("G", "G", 11, 24)] = 5006
    correct_expectation[("G", "G", 24, 37)] = 3784
    correct_expectation[("G", "G", 37, 37)] = 2046
    correct_expectation[("G", "G", 37, 2)] = 800  # Unique to G/G and to correct

    incorrect_expectation[("G", "A", 11, 24)] = 20
    incorrect_expectation[("G", "A", 24, 37)] = 20
    incorrect_expectation[("G", "A", 37, 37)] = 20
    incorrect_expectation[("G", "A", 2, 11)] = 10  # Unique to G/A and to incorrect
    incorrect_expectation[("G", "G", 11, 24)] = 20
    incorrect_expectation[("G", "G", 24, 37)] = 20
    incorrect_expectation[("G", "G", 37, 37)] = 20
    incorrect_expectation[("G", "G", 2, 11)] = 10

    # Run code and test
    correct_result, incorrect_result = _aggregate_mod_counts(
        correct, incorrect, rule="6bp"
    )

    assert correct_result == correct_expectation
    assert incorrect_result == incorrect_expectation


def test_error_prob_to_phred():

    # Case 1: Non-zero error probabilities
    array1 = np.array([0.001, 0.01, 0.1])
    expected_phreds1 = np.array([30, 20, 10])

    outcome1 = _error_prob_to_phred(array1)
    np.testing.assert_array_equal(outcome1, expected_phreds1)

    # Case 2: Some error probabilities are zero
    array2 = np.array([0.001, 0.0, 0.1, 0])
    expected_phreds2 = np.array([30, 60, 10, 60])

    outcome2 = _error_prob_to_phred(array2)
    np.testing.assert_array_equal(outcome2, expected_phreds2)

    # Case 3: All error probabilites are zero
    array3 = np.array([0.0, 0, 0.000000])
    expected_phreds3 = np.array([60, 60, 60])

    outcome3 = _error_prob_to_phred(array3)
    np.testing.assert_array_equal(outcome3, expected_phreds3)
