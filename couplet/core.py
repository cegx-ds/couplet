#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from couplet.align import (
    count_mismatches,
    needleman_wunsch,
    compute_mismatch_stats,
    get_aligned_record,
)
from couplet.resolve import resolve_reads
from couplet.trim import (
    trim_read,
    get_left_trimming_point,
    get_right_trimming_point,
    get_dynamic_right_trimming_point,
)
from enum import Enum

# Enum class describing the ways in which a resolved read may have been processed.
# The Enum values are used downstream in the form of custom .fastq tags.
class ResolutionTag(Enum):
    ACCEPTABLE = 1
    RESCUED = 2


def _get_record(original, alignment, missing_qual=0):
    """Expands a record according to an alignment string

    This function expands a record to include 'N' at positions in
    which there is a gap in the alignment string. For example,

    original  = (ATGC, FFFF)
    alignment = ATG-C
    new record = (ATGNC, FFF!F)
    """
    seq, phred = get_aligned_record(original, alignment, missing_qual=missing_qual)
    return SeqRecord(
        seq=Seq(seq),
        id=f"{original.id}",
        letter_annotations={"phred_quality": phred},
        description="",
    )


def _is_resolved(number_mismatches, L, mismatch_threshold):
    return True if number_mismatches <= mismatch_threshold * L else False


def resolve_read_pair(
    r1,
    r2,
    rule,
    resolve_phred_fn,
    mismatch_threshold,
    use_mismatch_aware_trimming=True,
    min_read_length=50,
    episeq_in_qname=False,
    log_additional_stats=True,
    orig_quals_in_sam_tag=False,
    xe_tag=False,
):
    L = min(len(r1), len(r2))

    # Count mismatches under naive alignment
    number_mismatches = count_mismatches(r1, r2, rule.permitted_pairs, match_n=True)
    stats = {"num_mismatch_naive_alignment": number_mismatches, "num_reads": 1}
    stats["input_bases"] = L

    # Case where reads r1 and r2 have sufficiently few mismatches that they can be naively aligned.
    if _is_resolved(number_mismatches, L, mismatch_threshold):
        result = resolve_reads(
            r1,
            r2,
            rule.rule_dict,
            rule.modifications,
            resolve_phreds=resolve_phred_fn,
            episeq_in_qname=episeq_in_qname,
            resolution_tag=ResolutionTag["ACCEPTABLE"].value,
            orig_quals_in_sam_tag=orig_quals_in_sam_tag,
            xe_tag=xe_tag,
        )
        stats["acceptable"] = 1
        stats[f"acceptable_length_{L}"] = 1
        stats["produced_bases"] = L

    # Case where r1 and r2 have too many mismatches and require  pairwise alignment.
    else:

        # Pairwise align reads and count mismatches
        align1, align2 = needleman_wunsch(r1.seq, r2.seq, rule)
        L = min(len(align2), len(align1))

        # Get mismatch positions. Include N's (e.g. from gaps) as mismatches.
        NW_mismatches = count_mismatches(
            align1, align2, rule.permitted_pairs, match_n=False, report_pos=True
        )

        # Perform trimming of N's. As a minimum we remove trailing N's using
        # the functions "left_trim" and "right_trim". Optionally, N's can be
        # trimmed from the right using a dynamic trimming process.
        left_trim = get_left_trimming_point(NW_mismatches)
        right_trim = get_right_trimming_point(NW_mismatches, L)

        if use_mismatch_aware_trimming:

            # If no dynamic trimming solution is found, right_trim is left unchanged
            NW_mismatches_right_trim = [m for m in NW_mismatches if m < right_trim]
            right_trim = get_dynamic_right_trimming_point(
                mismatches=NW_mismatches_right_trim,
                read_length=right_trim,
                min_read_length=min_read_length,
                left_trim=left_trim,
                threshold=mismatch_threshold,
            )

        # Get aligned record
        seq1, phred1 = get_aligned_record(r1, align1, missing_qual=0)
        seq2, phred2 = get_aligned_record(r2, align2, missing_qual=0)

        # Perform trimming and create SeqRecords
        result_rec1 = trim_read(r1.id, seq1, phred1, left_trim, right_trim)
        result_rec2 = trim_read(r2.id, seq2, phred2, left_trim, right_trim)

        # Recompute number of mismatches and L after trimming
        NW_mismatches = count_mismatches(
            result_rec1,
            result_rec2,
            rule.permitted_pairs,
            match_n=False,
            report_pos=True,
        )
        L = min(len(result_rec1), len(result_rec2))

        stats["num_mismatch_post_attempted_rescue"] = len(NW_mismatches)

        # Check if read is rescued after NW and trimming
        if len(NW_mismatches) <= L * mismatch_threshold and L >= min_read_length:

            result = resolve_reads(
                result_rec1,
                result_rec2,
                rule.rule_dict,
                rule.modifications,
                resolve_phreds=resolve_phred_fn,
                episeq_in_qname=episeq_in_qname,
                resolution_tag=ResolutionTag["RESCUED"].value,
                orig_quals_in_sam_tag=orig_quals_in_sam_tag,
                xe_tag=xe_tag,
            )
            stats["rescued"] = 1
            stats[f"rescued_length_{L}"] = 1
            stats["produced_bases"] = L

        else:
            result = None
            stats["discarded"] = 1
            stats[f"discarded_length_{L}"] = 1

    # Compute additional stats if requested
    if log_additional_stats:
        naive_mismatch_stats = compute_mismatch_stats(r1, r2, rule)

        if "acceptable" in stats:
            prefix = "acceptable"
        elif "rescued" in stats:
            prefix = "rescued"
        else:
            prefix = "discarded"

        # Add prefix to mismatch dictionary and merge into stats dictionary
        prefix_naive_mismatch_stats = {
            prefix + "_naive_" + key: value
            for key, value in naive_mismatch_stats.items()
        }
        stats = stats | prefix_naive_mismatch_stats

        # For all but acceptable reads, also record mismatch statistics post alignment
        if prefix in ["rescued", "discarded"]:
            aligned_mismatch_stats = compute_mismatch_stats(align1, align2, rule)

            # Add prefix to mismatch dictionary and merge into stats dictionary
            prefix_aligned_mismatch_stats = {
                prefix + "_aligned_" + key: value
                for key, value in aligned_mismatch_stats.items()
            }
            stats = stats | prefix_aligned_mismatch_stats

    return result, stats


def update_stats(stats, new_values):
    for k, v in new_values.items():
        if k not in stats:
            stats[k] = v
        else:
            stats[k] += v
    return stats
