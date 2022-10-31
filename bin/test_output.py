#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
from Bio import SeqIO
import gzip
import sys

fq1 = sys.argv[1]
fq2 = sys.argv[2]

for i, (record1, record2) in enumerate(
    zip(
        SeqIO.parse(gzip.open(fq1, "rt"), "fastq"),
        SeqIO.parse(gzip.open(fq2, "rt"), "fastq"),
    )
):
    assert (
        record1.id[:10] == record2.id[:10]
    ), f"Error on record {i}: {record1.id} <> {record2.id}"
    assert (
        record1.seq == record2.seq
    ), f"Error on record {i} ({record1.id}): {record1.seq} <> {record2.seq}"
    q1, q2 = (
        record1.letter_annotations["phred_quality"],
        record2.letter_annotations["phred_quality"],
    )
    assert q1 == q2, f"Error on record {i}: {q1} <> {q2}"
