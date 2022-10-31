#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
from Bio import SeqIO
from io import StringIO
from ssresolve.resolve import resolve_phred_prob
from ssresolve.rules import (
    five_bp_rule,
    five_bp_error_codes,
    five_bp_modifications,
    ResolutionRule,
)
from ssresolve.core import _is_resolved, resolve_read_pair


def _process_read_pair(
    rec1, rec2, fraction_allowed_mismatches, use_mismatch_aware_trimming
):
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

    result, stats = resolve_read_pair(
        rec1,
        rec2,
        rule=rule,
        resolve_phred_fn=resolve_phred_prob,
        mismatch_threshold=fraction_allowed_mismatches,
        min_read_length=1,
        use_mismatch_aware_trimming=use_mismatch_aware_trimming,
        log_additional_stats=True,
        xe_tag=True,
    )
    return result, stats


def test_needleman_wunsch_off_by_one():
    read_id = "A00536:97:HMTFKDMXX:1:1101:16767:1031"
    fq1 = "@" + read_id + "\n" "TNAATACTAA\n" "+\n" "F!FFFFFFFF"

    fq2 = "@" + read_id + "\n" "CTNAACACTA\n" "+\n" "FF!FFFFFFF"

    # R1 = TNAATACTAA
    # R2 = CTNAACACTA
    # A1 = -TNAATACTAA
    # A2 = CTNAATACTA-
    # GS = NTNAACACTAN
    # ES = NTNAACAPTAN
    #  Q = !B!BBBBBBB!
    rec1 = next(SeqIO.parse(StringIO(fq1), format="fastq"))
    rec2 = next(SeqIO.parse(StringIO(fq2), format="fastq"))
    result, _ = _process_read_pair(rec1, rec2, 1 / 9, False)
    assert result.seq == "TNAACACTA"
    assert result.id == read_id + " XE:Z:TNAACAPTA\tXR:i:2\tMM:Z:C+mh,1;"
    assert result.letter_annotations["phred_quality"] == [
        33,
        2,
        33,
        33,
        33,
        33,
        33,
        33,
        33,
    ]


def test_needleman_wunsch_and_mismatch_aware_trimming():

    # fmt: off
    read_id = "A00536:97:HMTFKDMXX:1:1101:16767:1031"
    fq1 = "@" + read_id +"\n" "AAAAACCCCGGCGGTTTTTC\n" "+\n" "FFFFFFFFFFF#FFFFFFF#"
    fq2 = "@" + read_id + "\n" "CAAAAACCCCAAAAATTATT\n" "+\n" "#FFFFFFFFFFFFFFFF#FF"
    # Resolution:                                     NAAAAAAAAAAANAAAANAAN
    # Indices:                                        012345678901234567890
    # fmt: on

    # Expectation:
    # - Pairwise align: Will translate sequences by one and add "N" to front of
    #   R1 and to back of R2.
    # - There will be mismatches at indices [0,12, 17, 20]
    # - Simple trimming will remove trailing mismatches --> [11,16]. This
    #   corresponds to a mismacth density of 2/19 = 0.105. If we require
    #   a density < 0.1, this would not pass. Thus we do dynamic trimming.
    # - Dynamic trimming will trim from right until mismatch at index 16
    # - There will then be one mismatch left in 16 bases. This is a density
    #   of 1/16 = 0.0625. If we require a density < 0.1, this should pass.

    rec1 = next(SeqIO.parse(StringIO(fq1), format="fastq"))
    rec2 = next(SeqIO.parse(StringIO(fq2), format="fastq"))
    result, _ = _process_read_pair(rec1, rec2, 0.1, True)

    assert result.seq == "AAAAACCCCGGNGGTT"

    assert result.letter_annotations["phred_quality"] == [33] * 11 + [2] + [33] * 4
    assert result.id == read_id + " XE:Z:AAAAAPPPPGG3GGTT\tXR:i:2\tMM:Z:C+mh,0,0,0,0;"


def test_is_resolved():
    number_mismatches = 1
    L = 100
    mismatch_threshold = 0.05
    status = _is_resolved(number_mismatches, L, mismatch_threshold)
    assert status

    number_mismatches = 8
    L = 100
    mismatch_threshold = 0.05
    status = _is_resolved(number_mismatches, L, mismatch_threshold)
    assert not status
