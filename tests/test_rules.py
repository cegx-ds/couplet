#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import numpy as np
import pytest

from ssresolve.rules import (
    get_permitted_pairs,
    create_nw_score_table,
    five_bp_rule,
    CHARLIST,
)


@pytest.fixture
def base_rule():
    # five_bp rule
    return {
        ("A", "A"): ("A", "A"),
        ("A", "C"): ("N", "0"),
        ("A", "G"): ("N", "1"),
        ("A", "T"): ("N", "2"),
        ("C", "A"): ("N", "3"),
        ("C", "C"): ("C", "P"),
        ("C", "G"): ("N", "4"),
        ("C", "T"): ("N", "5"),
        ("G", "A"): ("G", "G"),
        ("G", "C"): ("N", "6"),
        ("G", "G"): ("N", "n"),
        ("G", "T"): ("N", "7"),
        ("T", "A"): ("N", "8"),
        ("T", "C"): ("C", "C"),
        ("T", "G"): ("N", "9"),
        ("T", "T"): ("T", "T"),
    }


def test_get_permitted_pairs(base_rule):
    assert get_permitted_pairs(base_rule) == {
        ("T", "T"),
        ("A", "A"),
        ("T", "C"),
        ("C", "C"),
        ("G", "A"),
    }


def test_create_nw_score_table():
    score_table = create_nw_score_table(
        get_permitted_pairs(five_bp_rule),
        char_list=CHARLIST,
        match_award=1,
        gap_penalty=-2,
        mismatch_penalty=-1,
    )

    # Note that we assume (N,-) is a gap not a match
    expected_table = np.array(
        [
            # A   C   G   T   N   -
            [1, -1, -1, -1, 1, -2],  # A
            [-1, 1, -1, -1, 1, -2],  # C
            [1, -1, -1, -1, 1, -2],  # G
            [-1, 1, -1, 1, 1, -2],  # T
            [1, 1, 1, 1, 1, -2],  # N
            [-2, -2, -2, -2, -2, -2],  # N
        ]
    )

    np.testing.assert_array_equal(score_table, expected_table)
