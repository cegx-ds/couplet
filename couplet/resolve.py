#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import math
import functools
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import pysam

ILLEGAL_BP_PHRED_SCORE = 2


def resolve_bases(seq1, seq2, base_rule):
    """Resolve the bases in a pair of aligned sequences

    This functions assumes that a read pair has been aligned and any
    necessary preprocessing has been applied to the two reads.

    During resolution each pair of bases from the aligned positions across the
    two reads is passed to a rule to determine the resolved genetic and
    epigenetic symbols. Pairings absent from the rule are resolved to ('N', 'N')

    Args:
        seq1 (str): the first sequence
        seq2 (str): the second sequence
        base_rule (dict): a dictionary mapping pairs of bases to the resolved
            (genetic, epigenetic) pairs
    Returns:
        a str tuple (genomic_sequence, epigenomic_sequence) of sequences resolved
        using the rule
    """
    resolved = [base_rule.get((b1, b2), ("N", "N")) for b1, b2 in zip(seq1, seq2)]
    return ("".join([r[0] for r in resolved]), "".join([r[1] for r in resolved]))


def resolve_phred_min(qseq1, qseq2, seq1, seq2):
    """Resolve the quality scores in a pair of aligned quality sequences

    This functions assumes that a read pair has been aligned and any
    necessary preprocessing has been applied to the two reads.

    This function builds a new quality sequence made up of the lesser of the
    two quality scores associated with each base

    Args:
        qseq1 (str): the first read quality sequence
        qseq2 (str): the second read quality sequence


    Returns:
        an array representing the new quality sequence
    """

    return [min(q1, q2) for q1, q2 in zip(qseq1, qseq2)]


def resolve_phred_with_qtable(qseq1, qseq2, seq1, seq2, quality_table=None):
    """Resolve the quality scores in a pair of aligned quality sequences using a quality table

    This functions assumes that a read pair has been aligned and any
    necessary preprocessing has been applied to the two reads.

    This function builds a new quality sequence based on an input quality table.

    Args:
        qseq1 (str): the first read quality sequence
        qseq2 (str): the second read quality sequence
        r1 (str): the first read
        r2 (str): the second read
        quality_table: a phred quality table

    Returns:
        an array representing the new quality sequence
    """
    return [
        quality_table.get((b1, b2, q1, q2), ILLEGAL_BP_PHRED_SCORE)
        for b1, b2, q1, q2 in zip(seq1, seq2, qseq1, qseq2)
    ]


@functools.lru_cache(maxsize=1764)
def get_letter_joint_prob(Q1, Q2):
    """Return a new int based quality  on joint probability
    Where Q1 and Q2 are quality integers
    """

    R1_prob = 10.0 ** (-Q1 / 10.0)
    R2_prob = 10.0 ** (-Q2 / 10.0)

    prob = (R1_prob + R2_prob) - (R1_prob * R2_prob)

    jointQ = -10.0 * math.log10(prob)

    jointQ_int = int(
        jointQ
    )  # we are not adding 33 as biopython will add  https://github.com/biopython/biopython/blob/master/Bio/SeqIO/QualityIO.py
    retQ = min(42, jointQ_int)  # we dont want more better than K = 42 upper quality

    return retQ


def resolve_phred_prob(qseq1, qseq2, seq1, seq2):
    """Resolve the quality scores in a pair of aligned quality sequences

    This functions assumes that a read pair has been aligned and any
    necessary preprocessing has been applied to the two reads.

    This function builds a new quality sequence made up of the joint
    base call based on prob mode  the
    two quality scores associated with each base

    Args:
        qseq1 (str): the first read quality sequence
        qseq2 (str): the second read quality sequence

    Returns:
        an array representing the new quality sequence
    """

    return [get_letter_joint_prob(q1, q2) for q1, q2 in zip(qseq1, qseq2)]


def trim_n(genomic_seq, epigenomic_seq, resolved_phred, Q1, Q2):
    """Trim N's from the head and tail of a resolved read

    This function will remove N's from the head and the tail of a
    resolved genomic sequence and will remove the corresponding bases
    from the associated epigenomic sequence and the associated quality sequence

    Args:
        genomic_seq (str): the resolved genomic read
        epigenomic_seq (str): the resolved epigenomic read
        epigenomic_phred (str): the resolved phred sequence
        r_Q1 (list of int): the original Phred scores for the first read
        r_Q2 (list of int): the original Phred scores for the second read

    Returns:
        a str tuple (trimmed_genomic_sequence, trimmed_epigenomic_sequence,
                     trimmed_phred), trimmed_Q1, trimmed_Q2
    """
    # Strip N's from left
    l_genomic_seq = genomic_seq.lstrip("N")
    l_start = len(genomic_seq) - len(l_genomic_seq)
    l_epigenomic_seq = epigenomic_seq[l_start::]
    l_resolved_phred = resolved_phred[l_start::]
    l_Q1 = Q1[l_start::]
    l_Q2 = Q2[l_start::]

    # Strip N's from the right
    r_genomic_seq = l_genomic_seq.rstrip("N")
    r_end = len(l_genomic_seq) - len(r_genomic_seq)
    if r_end > 0:
        r_epigenomic_seq = l_epigenomic_seq[:-r_end:]
        r_resolved_phred = l_resolved_phred[:-r_end:]
        r_Q1 = l_Q1[:-r_end:]
        r_Q2 = l_Q2[:-r_end:]
    else:
        r_epigenomic_seq = l_epigenomic_seq
        r_resolved_phred = l_resolved_phred
        r_Q1 = l_Q1
        r_Q2 = l_Q2

    return (r_genomic_seq, r_epigenomic_seq, r_resolved_phred, r_Q1, r_Q2)


def get_fragment_length_tag(description):

    """Takes 'description' field from a SeqRecord object and checks to 
    see if it contains the fragment length tag 'xl:...'. 
    If so, returns the tag.

    Args:
        description (str): the 'description' field of a SeqRecord
    Returns:
        fragment_length_tag (str): "" if none present
    """

    xl_tag = re.search(r'xl:i:\d*',description)
    
    if xl_tag is not None:
        xl_tag = xl_tag.group()
        fragment_length_tag = f"\t{xl_tag}"
    else:
        fragment_length_tag = ""

    return fragment_length_tag


def resolve_reads(
    r1,
    r2,
    base_rule,
    modifications,
    resolve_phreds,
    episeq_in_qname=False,
    resolution_tag=None,
    orig_quals_in_sam_tag=False,
    xe_tag=False,
):
    """Resolve a pair of aligned reads

    This functions assumes that a read pair has been aligned and any necessary
    preprocessing has been applied to the two reads.

    The following steps happen during resolution:
    1. Each pair of bases from the aligned positions across the two reads is passed to
       a rule to determine the resolved genetic and epigenetic symbols
    2. Each pair of phred scores from the aligned positions across the two reads is
       processed to determine the resolved phred score. This is currently limited to the
       most basic resolution logic, where we take the minimum of the two phred scores

    Args:
        r1 (SeqIO.SeqRecord)         : the first read
        r2 (SeqIO.SeqRecord)         : the second read
        base_rule (dict)             : a dictionary mapping pairs of bases to the
                                       resolved (genetic, epigenetic) pairs
        resolve_phreds (function)    : a function for resolving phred scores
        episeq_in_qname              : A boolean denoting if the epigenetic sequence
                                       should be stored inside the fastq read name.
                                       Default is False.
        resolution_tag               : An integer describing the processing method (e.g.
                                       acceptable or rescued reads) which is to be added
                                       to the custom information for the resolved read.
                                       The tag is included as a SAM flag with the symbol
                                       "R" (for resolution). E.g. "XR:i:1". If this
                                       input is None, the tag won't be added.
        orig_quals_in_sam_tag (bool) : A boolean determining whether the original
                                       qualities (Q1, Q2) should be stored in a SAM
                                       tag.
        xe_tag (bool)                : A boolean determining whether base modification
                                       information should be stored in a custom XE tag
    Returns:
        a SeqIO.SeqRecord corresponding to the resolved read
    """

    (genomic_seq, epigenomic_seq) = resolve_bases(r1.seq, r2.seq, base_rule)

    # Get quality scores for R1 and R2
    Q1 = r1.letter_annotations["phred_quality"]
    Q2 = r2.letter_annotations["phred_quality"]

    # If len(R1) != len(R2), the resolved sequence will be shorter than the
    # Q-scores. Thus we trim the Q-scores to match the resolved sequence.
    Q1 = Q1[: len(genomic_seq)]
    Q2 = Q2[: len(genomic_seq)]

    resolved_phred = resolve_phreds(
        Q1,
        Q2,
        r1.seq,
        r2.seq,
    )

    # TODO: Fix this
    for i, c in enumerate(genomic_seq):
        if c == "N":
            resolved_phred[i] = ILLEGAL_BP_PHRED_SCORE

    # TODO: This is only relevant for naively aligned reads.
    # TODO: This should be handled in core.py, not here.
    # Trim N's from either end
    (genomic_seq, epigenomic_seq, resolved_phred, Q1, Q2) = trim_n(
        genomic_seq, epigenomic_seq, resolved_phred, Q1, Q2
    )

    # Define how to store the epigenetic sequence. Either in read name or as separate field.
    if episeq_in_qname:
        episeq_formatted = f":cseq:{epigenomic_seq}"
    else:
        episeq_formatted = f" XE:Z:{epigenomic_seq}"

    # Add a resolution tag describing how the read was produced
    if resolution_tag is None:
        resolution_tag = ""
    else:
        resolution_tag = f"\tXR:i:{resolution_tag}"

    # Add a base modification tag
    base_mod_tag = compute_base_mod_tag(genomic_seq, epigenomic_seq, modifications)

    # Export original quality scores if requested
    if orig_quals_in_sam_tag is True:

        # Convert Phred scores from integer encoding (2,30) to symbolic encoding (#,?)
        Q1_string = pysam.qualities_to_qualitystring(Q1)
        Q2_string = pysam.qualities_to_qualitystring(Q2)
        orig_quals_tag = f"\tXQ:Z:{Q1_string}\tYQ:Z:{Q2_string}"

    else:
        orig_quals_tag = ""

    if not xe_tag:
        episeq_formatted = ""

    fragment_length_tag = get_fragment_length_tag(r1.description)

    result = SeqRecord(
        seq=Seq(genomic_seq),
        id=f"{r1.id}{episeq_formatted}{resolution_tag}{base_mod_tag}"
        f"{orig_quals_tag}{fragment_length_tag}",
        letter_annotations={"phred_quality": resolved_phred},
        description="",
    )

    return result


def compute_base_mod_tag(genomic_seq, epigenomic_seq, modifications):

    """Computes a SAM tag for base modifications.

    This function takes genomic and epigenomic sequences and computes a SAM tag
    associated with a set of input modifications. Modifications are assumed to be of
    the form: { "epigenomic code": [ "genomic base", "SAM code" ] }. The epigenomic
    code can be a single base (e.g. "P") or a pair of bases (e.g. "PQ" or "PG"). The
    genomic base must be a single base, representing the genomic base upon which the
    modification acts. The SAM code must be one of the modification codes given in the
    SAM specificiation.

    The code returns a tag of the form: "\tMM:Z:C+m,5,3,0,20;C+h,10,14,3;", i.e.
    each modification is separated by a semicolon.


    Args:
        genomic_seq    : A sequence of genomic bases (A, C, G, T)
        epigenomic_seq : A sequence made up of genomic bases (A, C, G, T), epigenomic
                         bases (e.g. P and Q) and error codes (0, 1, 2, ...).
        modifications  : A dictionary of modifications. These are of the form:
                         { "epigenomic code": [ "genomic base", "SAM code" ] } e.g.
                         { "PQ": [ "C", "m" ], "PG": [ "C", "h"] }

    Returns:
        A tag with base modification information.
    """

    tag = "\tMM:Z:"

    # Loop over modifications
    for key, values in modifications.items():

        # Example: mod = { "P": ["C", "mh"] }
        epi_context = key
        gen_base = values[0]
        SAM_code = values[1]

        # Create initial part of tag for current modification
        tag += f"{gen_base}+{SAM_code}"

        # Counts number of times since a modification has been seen
        non_mod_counter = 0

        # Order of terms separated by OR in the regex defines preference, e.g.
        # re.findall("C|PG|P") means that we find "PG" before we find "P".
        # If { "P": ["C", "mh"] }, then find epi_context and epi_context[0] still return
        # all relevant results and the 1-let correct match triggers tag modification
        for match in re.findall(
            f"{gen_base}|{epi_context}|{epi_context[0]}", epigenomic_seq
        ):
            if match == epi_context:
                tag += "," + str(non_mod_counter)
                non_mod_counter = 0
            else:
                non_mod_counter += 1

        tag += ";"

    return tag
