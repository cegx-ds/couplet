#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import pytest

from io import StringIO
from Bio import SeqIO
import numpy as np
from numba.typed import List

from ssresolve.align import count_mismatches, compute_mismatch_stats, needleman_wunsch
from ssresolve.align import _nw_create_matrix_numba as nw_create_matrix_numba
from ssresolve.rules import (
    five_bp_rule,
    five_bp_error_codes,
    five_bp_modifications,
    ResolutionRule,
)


@pytest.fixture
def permitted_pairs():
    return {("A", "A"), ("C", "C"), ("G", "T")}


def test_align(permitted_pairs):
    r1 = "ACTG"
    r2 = "ACTT"
    assert count_mismatches(r1, r2, permitted_pairs) == 1
    assert count_mismatches(r1, r2, permitted_pairs, report_pos=True) == [2]

    r1 = "ACG"
    r2 = "ACT"
    assert count_mismatches(r1, r2, permitted_pairs) == 0
    assert count_mismatches(r1, r2, permitted_pairs, report_pos=True) == []


def test_different_orders(permitted_pairs):
    r1 = "ACTGT"
    r2 = "ACTTG"
    assert count_mismatches(r1, r2, permitted_pairs) == 2
    assert count_mismatches(r1, r2, permitted_pairs, report_pos=True) == [2, 4]


def test_different_lengths(permitted_pairs):
    r1 = "ACTGT"
    r2 = "ACTTGT"
    assert count_mismatches(r1, r2, permitted_pairs) == 2
    assert count_mismatches(r1, r2, permitted_pairs, report_pos=True) == [2, 4]


def test_match_n(permitted_pairs):
    r1 = "ACTGT"
    r2 = "ACTTN"
    assert count_mismatches(r1, r2, permitted_pairs, match_n=True) == 1
    assert count_mismatches(r1, r2, permitted_pairs, match_n=False) == 2
    assert count_mismatches(r1, r2, permitted_pairs, match_n=True, report_pos=True) == [
        2
    ]
    assert count_mismatches(
        r1, r2, permitted_pairs, match_n=False, report_pos=True
    ) == [2, 4]


@pytest.fixture
def rule():
    gap_penalty = -2
    match_award = 1
    mismatch_penalty = -1
    end_gap_penalty = -1
    rule = ResolutionRule(
        five_bp_rule,
        five_bp_error_codes,
        five_bp_modifications,
        match_award=match_award,
        gap_penalty=gap_penalty,
        mismatch_penalty=mismatch_penalty,
        end_gap_penalty=end_gap_penalty,
    )
    return rule


def test_needleman_wunsch_identical(rule):
    assert needleman_wunsch("TGAGTTT", "TGAGTTT", rule) == ("TGAGTTT", "TGAGTTT")


def test_needleman_wunsch_substitution(rule):
    # Seqs unchanged if substitution only
    assert needleman_wunsch("TCGA", "TCCA", rule) == ("TCGA", "TCCA")


def test_needleman_wunsch_shift(rule):
    # Second read is shifted by 1
    assert needleman_wunsch("TCCATCCATCCA", "CCATCCATCCAT", rule) == (
        "TCCATCCATCCA-",
        "-CCATCCATCCAT",
    )


def test_needleman_wunsch_gap_at_end(rule):
    # Prefer to have gaps at the end when gap_penalty < end_gap_penalty
    assert needleman_wunsch("TCAATCAATCA", "CAATCAATCAA", rule) == (
        "TCAATCAATCA-",
        "-CAATCAATCAA",
    )
    assert needleman_wunsch("GAATTT", "AAAC", rule) == (
        "GAATTT",
        "AAAC--",
    )

    rule.end_gap_penalty = -2
    # When gap_penalty = end_gap_penalty, we expect matches to be
    # consumed before gaps when traversing the matrix from bottom-right
    assert needleman_wunsch("GAATTT", "AAAC", rule) == (
        "GAATTT",
        "AAA--C",
    )

    rule.end_gap_penalty = -5
    # Prefer not to have gaps at the end when gap_penalty > end_gap_penalty
    assert needleman_wunsch("TCAATCAATCA", "CAATCAATCAA", rule) == (
        "TCAATCAATC-A",
        "C-AATCAATCAA",
    )
    assert needleman_wunsch("GAATTT", "AAAC", rule) == (
        "GAATTT",
        "AAA--C",
    )


def test_score_matrix(rule):
    s1 = List([rule.index[c] for c in "GAATTT"])
    s2 = List([rule.index[c] for c in "AAAC"])
    score = nw_create_matrix_numba(
        s1,
        s2,
        rule.nw_score_table,
        gap_penalty=rule.gap_penalty,
        end_gap_penalty=rule.end_gap_penalty,
    )

    expected_score = np.array(
        [
            [0.0, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0],
            [-1.0, 1.0, 0.0, -1.0, -3.0, -5.0, -6.0],
            [-2.0, 0.0, 2.0, 1.0, -1.0, -3.0, -4.0],
            [-3.0, -1.0, 1.0, 3.0, 1.0, -1.0, -2.0],
            [-4.0, -2.0, 0.0, 2.0, 4.0, 3.0, 2.0],
        ]
    )

    np.testing.assert_array_equal(score, expected_score)


def test_compute_mismatch_stats(rule):
    r1 = "ACTGTAAGGA"
    r2 = "ACTTNCAGAC"
    expected = {
        "7_pos_3": 1,
        "0_pos_5": 1,
        "n_pos_7": 1,
        "0_pos_9": 1,
        "7": 1,
        "0": 2,
        "n": 1,
    }
    assert compute_mismatch_stats(r1, r2, rule) == expected
