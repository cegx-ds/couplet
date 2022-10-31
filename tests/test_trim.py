#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import pytest
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from io import StringIO

from ssresolve.trim import (
    get_dynamic_right_trimming_point,
    trim_read,
    get_left_trimming_point,
    get_right_trimming_point,
)
from ssresolve.align import needleman_wunsch


def test_trim_read():

    # Trim read from 16 bases down to 10
    read_id = "this:is:a:read:id"
    seq = "ACGT" * 4
    phred = [31, 32, 33, 34] * 4
    trim1 = 0
    trim2 = 10

    outcome_rec = trim_read(read_id, seq, phred, trim1, trim2)

    assert outcome_rec.seq == "ACGTACGTAC"
    assert outcome_rec.id == read_id
    assert outcome_rec.letter_annotations["phred_quality"] == [
        31,
        32,
        33,
        34,
        31,
        32,
        33,
        34,
        31,
        32,
    ]
    assert outcome_rec.description == ""


def test_get_left_trimming_point():

    # Test no mismatches provided
    assert get_left_trimming_point([]) == 0

    # Test no need for trimming
    assert get_left_trimming_point([1]) == 0
    assert get_left_trimming_point([80, 90, 100]) == 0

    # Test trimming works
    assert get_left_trimming_point([0, 1, 2, 3, 10]) == 4

    # Test all mismatches are left trimmed
    assert get_left_trimming_point([0, 1, 2, 3, 4]) == 5


def test_get_right_trimming_point():

    # Test no mismatches provided
    read_length = 10
    assert get_right_trimming_point([], read_length) == read_length

    # Test no need for trimming
    assert get_right_trimming_point([0], read_length) == read_length
    assert get_right_trimming_point([5, 7], read_length) == read_length

    # Test trimming works
    assert get_right_trimming_point([7, 8, 9], read_length) == 7
    assert get_right_trimming_point([9], read_length) == 9

    # Test all bases are mismatches
    assert get_right_trimming_point([0, 1, 2, 3, 4], 5) == 0


def test_get_dynamic_right_trimming_point():

    # Test no mismatches provided
    assert get_dynamic_right_trimming_point([], 100, 10, 0) == 100

    # Test no trimming required
    assert get_dynamic_right_trimming_point([20, 40, 50], 100, 10, 0) == 100

    # Test dynamic trimming
    assert (
        get_dynamic_right_trimming_point(
            [20, 40, 50, 70, 80, 85, 90], 100, 10, 0, threshold=0.05
        )
        == 80
    )
    assert (
        get_dynamic_right_trimming_point(
            [0, 1, 2, 3, 4, 70, 80, 85, 90, 95, 99], 100, 10, 5, threshold=0.05
        )
        == 95
    )

    # Test a case where all the N's are left trailing
    assert (
        get_dynamic_right_trimming_point(
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 100, 10, 11, threshold=0.05
        )
        == 100
    )

    # Test that minimum read length is used
    assert (
        get_dynamic_right_trimming_point(
            [10, 20, 30, 40, 50, 60, 70, 80, 90], 100, 30, 0, threshold=0.05
        )
        == 100
    )

    # Test a case where read_length is zero (which can arise when all bases are mismatches)
    assert (
        get_dynamic_right_trimming_point(
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 0, 5, 11, threshold=0.05
        )
        == 0
    )
