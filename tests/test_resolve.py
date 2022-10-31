#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
from collections import namedtuple
import io

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest
import functools


from ssresolve.resolve import resolve_reads, get_letter_joint_prob
from ssresolve.resolve import resolve_phred_with_qtable
from ssresolve.resolve import resolve_phred_min
from ssresolve.resolve import compute_base_mod_tag
from ssresolve.rules import five_bp_modifications
from ssresolve.rules import six_bp_modifications


DummyRecord = namedtuple("DummyRecord", "id seq letter_annotations")

# Subset of five_bp rule
@pytest.fixture
def rule():
    return {("A", "A"): ("A", "A"), ("C", "C"): ("C", "P"), ("G", "A"): ("G", "G")}


def test_resolve_reads(rule):

    # Case 1: N at end of read
    r1 = SeqRecord(
        seq=Seq("CGN"), id="@ABC1", letter_annotations={"phred_quality": [31, 29, 0]}
    )
    r2 = SeqRecord(
        seq=Seq("CAA"), id="@ABC2", letter_annotations={"phred_quality": [30, 32, 33]}
    )

    result = resolve_reads(
        r1, r2, rule, five_bp_modifications, resolve_phred_min, xe_tag=False
    )
    result_episeq_in_qname = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        episeq_in_qname=True,
        xe_tag=False,
    )
    result_resolution_tag = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        resolution_tag=1,
        xe_tag=False,
    )

    result_orig_quals_in_sam_tag = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        orig_quals_in_sam_tag=True,
        xe_tag=True,
    )
    result_xe_tag = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        xe_tag=True,
    )

    result_all_tags = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        resolution_tag=2,
        orig_quals_in_sam_tag=True,
        xe_tag=True,
    )

    assert result.seq == "CG"
    assert result.letter_annotations["phred_quality"] == [30, 29]
    assert result.id == "@ABC1\tMM:Z:C+mh,0;"
    assert result_episeq_in_qname.id == "@ABC1\tMM:Z:C+mh,0;"
    assert result_resolution_tag.id == "@ABC1\tXR:i:1\tMM:Z:C+mh,0;"
    assert (
        result_orig_quals_in_sam_tag.id
        == "@ABC1 XE:Z:PG\tMM:Z:C+mh,0;\tXQ:Z:@>\tYQ:Z:?A"
    )
    assert result_xe_tag.id == "@ABC1 XE:Z:PG\tMM:Z:C+mh,0;"
    assert result_all_tags.id == "@ABC1 XE:Z:PG\tXR:i:2\tMM:Z:C+mh,0;\tXQ:Z:@>\tYQ:Z:?A"

    # Case 2: Read 1 is longer than read 2
    r1 = SeqRecord(
        seq=Seq("CGAAAA"),
        id="@ABC1",
        letter_annotations={"phred_quality": [31, 29, 30, 30, 30, 30]},
    )
    r2 = SeqRecord(
        seq=Seq("CA"), id="@ABC2", letter_annotations={"phred_quality": [30, 32]}
    )

    result = resolve_reads(
        r1, r2, rule, five_bp_modifications, resolve_phred_min, xe_tag=False
    )
    result_episeq_in_qname = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        episeq_in_qname=True,
        xe_tag=False,
    )
    result_resolution_tag = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        resolution_tag=1,
        xe_tag=False,
    )

    result_orig_quals_in_sam_tag = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        orig_quals_in_sam_tag=True,
        xe_tag=False,
    )

    result_xe_tag = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        resolution_tag=1,
        xe_tag=True,
    )

    result_all_tags = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        resolution_tag=2,
        orig_quals_in_sam_tag=True,
        xe_tag=True,
    )

    assert result.seq == "CG"
    assert result.letter_annotations["phred_quality"] == [30, 29]
    assert result.id == "@ABC1\tMM:Z:C+mh,0;"
    assert result_episeq_in_qname.id == "@ABC1\tMM:Z:C+mh,0;"
    assert result_resolution_tag.id == "@ABC1\tXR:i:1\tMM:Z:C+mh,0;"
    assert result_orig_quals_in_sam_tag.id == "@ABC1\tMM:Z:C+mh,0;\tXQ:Z:@>\tYQ:Z:?A"
    assert result_xe_tag.id == "@ABC1 XE:Z:PG\tXR:i:1\tMM:Z:C+mh,0;"
    assert result_all_tags.id == "@ABC1 XE:Z:PG\tXR:i:2\tMM:Z:C+mh,0;\tXQ:Z:@>\tYQ:Z:?A"


def test_resolve_reads_with_qtable(rule):

    # Expectations
    # - First two bases will be resolved according to the Phred table
    # - Third base does not exist in the table, so will be given Phred of 2
    # - Fourth base is a trailing N, so will be trimmed away
    r1 = SeqRecord(
        seq=Seq("CGAN"),
        id="@ABC1",
        letter_annotations={"phred_quality": [31, 29, 12, 2]},
    )
    r2 = SeqRecord(
        seq=Seq("CAAA"),
        id="@ABC2",
        letter_annotations={"phred_quality": [30, 32, 33, 33]},
    )

    phred_table = {
        ("C", "C", 31, 30): 35,
        ("G", "A", 29, 32): 32,
        ("A", "A", 28, 33): 31,
    }

    resolve_phreds = functools.partial(
        resolve_phred_with_qtable, quality_table=phred_table
    )

    result = resolve_reads(
        r1, r2, rule, five_bp_modifications, resolve_phreds, xe_tag=False
    )

    result_episeq_in_qname = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phreds,
        episeq_in_qname=True,
        xe_tag=False,
    )
    result_resolution_tag = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phreds,
        resolution_tag=1,
        xe_tag=False,
    )

    assert result.seq == "CGA"
    assert result.letter_annotations["phred_quality"] == [35, 32, 2]
    assert result.id == "@ABC1\tMM:Z:C+mh,0;"
    assert result_episeq_in_qname.id == "@ABC1\tMM:Z:C+mh,0;"
    assert result_resolution_tag.id == "@ABC1\tXR:i:1\tMM:Z:C+mh,0;"


def test_trim_n(rule):
    r1 = SeqRecord(
        seq=Seq("NNNAAAAAAANN"),
        id="@ABC1",
        letter_annotations={
            "phred_quality": [0, 0, 0, 31, 29, 31, 29, 31, 29, 29, 0, 0]
        },
    )
    r2 = SeqRecord(
        seq=Seq("NCTAAAAAAANG"),
        id="@ABC2",
        letter_annotations={
            "phred_quality": [0, 31, 31, 31, 29, 31, 29, 31, 29, 29, 0, 31]
        },
    )

    result = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        orig_quals_in_sam_tag=True,
        xe_tag=False,
    )

    assert result.seq == "AAAAAAA"
    assert result.letter_annotations["phred_quality"] == [31, 29, 31, 29, 31, 29, 29]
    assert result.id == "@ABC1\tMM:Z:C+mh;\tXQ:Z:@>@>@>>\tYQ:Z:@>@>@>>"


def test_trim_n_no_n(rule):
    r1 = SeqRecord(
        seq=Seq("AAAAAAA"),
        id="@ABC1",
        letter_annotations={"phred_quality": [31, 29, 31, 29, 31, 29, 29]},
    )
    r2 = SeqRecord(
        seq=Seq("AAAAAAA"),
        id="@ABC2",
        letter_annotations={"phred_quality": [31, 29, 31, 29, 31, 29, 29]},
    )

    result = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        orig_quals_in_sam_tag=True,
        xe_tag=False,
    )

    assert result.seq == "AAAAAAA"
    assert result.letter_annotations["phred_quality"] == [31, 29, 31, 29, 31, 29, 29]
    assert result.id == "@ABC1\tMM:Z:C+mh;\tXQ:Z:@>@>@>>\tYQ:Z:@>@>@>>"


def test_write_read(rule):
    r1 = SeqRecord(
        seq=Seq("CGN"), id="ABC1", letter_annotations={"phred_quality": [31, 29, 0]}
    )
    r2 = SeqRecord(
        seq=Seq("CAA"), id="ABC2", letter_annotations={"phred_quality": [30, 32, 33]}
    )

    result = resolve_reads(
        r1,
        r2,
        rule,
        five_bp_modifications,
        resolve_phred_min,
        resolution_tag=1,
        orig_quals_in_sam_tag=True,
        xe_tag=True,
    )
    out_file = io.StringIO()
    SeqIO.write(result, out_file, "fastq")
    print(result)
    print(out_file.getvalue())
    assert (
        out_file.getvalue()
        == """@ABC1 XE:Z:PG\tXR:i:1\tMM:Z:C+mh,0;\tXQ:Z:@>\tYQ:Z:?A\nCG\n+\n?>\n"""
    )


def test_get_letter_joint_prob(Q1=37, Q2=37):
    assert 33 == get_letter_joint_prob(Q1, Q2)


def test_compute_base_mod_tag():

    # Case 1: Mix of C and modC in 5-Letter mode
    gen_seq1 = "ACTACGACCCCTCGACGAACA"
    epi_seq1 = "ACTAPGACCCCTPGAPGAACA"
    tag1 = compute_base_mod_tag(gen_seq1, epi_seq1, five_bp_modifications)
    expected_tag1 = "\tMM:Z:C+mh,1,4,0;"

    assert tag1 == expected_tag1

    # Case 2: No modifications in 5-Letter mode
    gen_seq2 = "ACAACAA"
    epi_seq2 = "ACAACAA"
    tag2 = compute_base_mod_tag(gen_seq2, epi_seq2, five_bp_modifications)
    expected_tag2 = "\tMM:Z:C+mh;"

    assert tag2 == expected_tag2

    # Case 3: All C's are modified in 5-Letter mode
    gen_seq3 = "CGCGCG"
    epi_seq3 = "PGPGPG"
    tag3 = compute_base_mod_tag(gen_seq3, epi_seq3, five_bp_modifications)
    expected_tag3 = "\tMM:Z:C+mh,0,0,0;"

    assert tag3 == expected_tag3

    # Case 4: Modifications at beginning and end in 5-Letter mode
    gen_seq4 = "CGAACCAACG"
    epi_seq4 = "PGAACCAAPG"
    tag4 = compute_base_mod_tag(gen_seq4, epi_seq4, five_bp_modifications)
    expected_tag4 = "\tMM:Z:C+mh,0,2;"

    assert tag4 == expected_tag4

    # Case 5: C at ends in 5-Letter mode
    gen_seq5 = "CAACGAAC"
    epi_seq5 = "CAAPGAAC"
    tag5 = compute_base_mod_tag(gen_seq5, epi_seq5, five_bp_modifications)
    expected_tag5 = "\tMM:Z:C+mh,1;"

    assert tag5 == expected_tag5

    # Case 6: Three types of modifications in 6-Letter mode
    gen_seq6 = "ACGTCAACGCQTCCGAAC"
    epi_seq6 = "ACGTPAAPGPQTCPGAAP"
    tag6 = compute_base_mod_tag(gen_seq6, epi_seq6, six_bp_modifications)
    expected_tag6 = "\tMM:Z:C+m,3;C+h,2,2;C+mh,1,0,0,1,0;"

    assert tag6 == expected_tag6

    # Case 7: C and hmC modifications (no mC) in 6-Letter mode
    gen_seq7 = "ACGTCAACGCQTCCGAAC"
    epi_seq7 = "ACGTPAAPGPGTCPGAAP"
    tag7 = compute_base_mod_tag(gen_seq7, epi_seq7, six_bp_modifications)
    expected_tag7 = "\tMM:Z:C+m;C+h,2,0,1;C+mh,1,0,0,1,0;"

    assert tag7 == expected_tag7
