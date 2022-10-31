#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import itertools as it
import numpy as np
from numba import jit
from numba.typed import List
from collections import Counter


def count_mismatches(r1, r2, permitted_pairs, match_n=True, report_pos=False):
    """Count the number of mismatched bases if we use the naive algorithm
    Reads are aligned naively by matching the i-th base in r1 with the i-th base
    in r2. This function is used to assess whether or not a pair of reads can be aligned
    by counting the number of base pairs which are not permitted under the selected
    rule. This function ignores overhanging bases in cases where one read is longer
    than the other.
    Args:
        r1 (list-like): the first read
        r2 (list-like): the second read
        permitted_pairs (iterable): the set of pairs which are permitted by the rule
        match_n (bool): whether we should treat N as an automatic match
        report_pos (bool): whether to report the position (indices) of mismatches
            instead of the number of mismatches
    Returns:
        The number of aligned base pairs which are not permitted by the rule
    """
    if match_n:
        matcher = (
            lambda a, b: (a, b) in permitted_pairs or "N" in [a, b] or "-" in [a, b]
        )
    else:
        matcher = lambda a, b: (a, b) in permitted_pairs

    # Either return the position of mismatches or the total count
    if report_pos:
        length = min(len(r1), len(r2))  # Only consider shared bases
        return [index for index in range(length) if not matcher(r1[index], r2[index])]
    else:
        return sum(1 for (a, b) in zip(r1, r2) if not matcher(a, b))


def compute_mismatch_stats(r1, r2, rule):

    """Returns a dictionary with error code positions and counts from resolving
    r1 and r2.

    This function takes two reads, r1 and r2, and extracts information about
    error codes resulting from pairwise resolution. Error codes and rules are
    defined in the ResolutionRule (rule) object.

    Args:
        r1 (iterable): The first read
        r2 (iterable): The second red
        rule: ResolutionRule object describing error codes and base rule

    Returns:
        A dictionary with the total count and positions of error codes arising
        from the resolution of the input reads. The output looks like:
        { "3_pos_89": 1, "3_pos_102": 1, "3": 2 }
    """
    error_codes = rule.error_codes
    base_rule = rule.rule_dict

    # E.g. A, 0, 5, P
    # Get epigenetic codes when resolving r1, r2
    # E.g. [ "A", "C", "C", "0", "T", "G", "A", "P", "3", "A", "8", "C" ]
    epi_codes = [base_rule.get((a, b), ("None", "None"))[1] for (a, b) in zip(r1, r2)]

    # Extract error codes and positions
    error_info = [
        (code, pos) for (pos, code) in enumerate(epi_codes) if code in error_codes
    ]

    # Create error code dictionaries with positional and count information, e.g.
    # { '0': 3, '1': 2, '3': 0, '4': 11, ...}
    # { '0_pos_30': 1, '0_pos_75': 1, '0_pos_89': 1, '1_pos_13': 1, ...}
    counts_dict = Counter([c[0] for c in error_info])
    pos_dict = {f"{code}_pos_{pos}": 1 for (code, pos) in error_info}

    # Dictionary merge operator (|) requires Python >= 3.9
    return counts_dict | pos_dict


@jit(nopython=True)
def _nw_create_matrix_numba(seq1, seq2, score_table, gap_penalty, end_gap_penalty):
    # Acknowledgement https://wilkelab.org/classes/SDS348/2019_spring/labs/lab13-solution.html

    # Store length of two sequences.
    # Note that n will be the number of columns and m the number of rows in the NW score matrix
    n, m = len(seq1), len(seq2)

    # Generate matrix of zeros to store scores
    score = np.zeros((m + 1, n + 1))

    # Calculate score table
    # Fill out first column and first row
    score[:, 0] = end_gap_penalty * np.arange(m + 1)
    score[0, :] = end_gap_penalty * np.arange(n + 1)

    # Fill out all other values in the score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if i == m or j == n:
                gap_penalty_fill = end_gap_penalty
            else:
                gap_penalty_fill = gap_penalty
            # Calculate the score by checking the top, left, and diagonal cells
            match = score[i - 1][j - 1] + score_table[(seq1[j - 1], seq2[i - 1])]
            delete = score[i - 1][j] + gap_penalty_fill
            insert = score[i][j - 1] + gap_penalty_fill

            # Record the maximum score from the three possible scores calculated above
            score[i][j] = max(match, delete, insert)

    return score


@jit(nopython=True)
def _nw_get_alignment_numba(
    seq1, seq2, score_table, score, gap_idx, gap_penalty, end_gap_penalty
):

    n, m = len(seq1), len(seq2)

    # Traceback and compute the alignment
    # Create variables to store alignment
    align1 = []
    align2 = []

    # Start from the bottom right cell in matrix
    i = m
    j = n

    # We'll use i and j to keep track of where we are in the matrix, just like above
    while i > 0 and j > 0:  # end touching the top or the left edge
        score_current = score[i][j]
        score_diagonal = score[i - 1][j - 1]
        score_up = score[i][j - 1]
        score_left = score[i - 1][j]

        if i == m or j == n:
            gap_penalty_fill = end_gap_penalty
        else:
            gap_penalty_fill = gap_penalty

        # Check to figure out which cell the current score was calculated from,
        # then update i and j to correspond to that cell.
        temp_mval = score_table[(seq1[j - 1], seq2[i - 1])]

        if score_current == score_diagonal + temp_mval:
            # total_score += score_current
            align1.append(seq1[j - 1])
            align2.append(seq2[i - 1])
            i -= 1
            j -= 1
        elif score_current == score_up + gap_penalty_fill:
            # total_score += score_current
            align1.append(seq1[j - 1])
            align2.append(gap_idx)
            j -= 1
        elif score_current == score_left + gap_penalty_fill:
            # total_score += score_current
            align1.append(gap_idx)
            align2.append(seq2[i - 1])
            i -= 1
        else:
            # We should never end up here, but if we did we would be stuck in
            # an infinite loop.
            assert (
                False
            ), "Unreachable statement in couplet.align._nw_get_alignment_numba."

    # Finish tracing up to the top left cell
    while j > 0:
        align1.append(seq1[j - 1])
        align2.append(gap_idx)
        j -= 1
    while i > 0:
        align1.append(gap_idx)
        align2.append(seq2[i - 1])
        i -= 1

    # Since we traversed the score matrix from the bottom right, our two sequences will be reversed.
    # These two lines reverse the order of the characters in each sequence.
    align1 = align1[::-1]
    align2 = align2[::-1]

    return align1, align2


def needleman_wunsch(seq1, seq2, rule):
    s1 = List([rule.index[c] for c in seq1])
    s2 = List([rule.index[c] for c in seq2])
    score = _nw_create_matrix_numba(
        s1,
        s2,
        rule.nw_score_table,
        gap_penalty=rule.gap_penalty,
        end_gap_penalty=rule.end_gap_penalty,
    )
    align1, align2 = _nw_get_alignment_numba(
        s1,
        s2,
        rule.nw_score_table,
        score,
        gap_idx=rule.char_list.index("-"),
        gap_penalty=rule.gap_penalty,
        end_gap_penalty=rule.end_gap_penalty,
    )
    align1 = "".join([rule.char_list[i] for i in align1])
    align2 = "".join([rule.char_list[i] for i in align2])
    return align1, align2


def get_aligned_record(record, alignment, missing_qual):
    original = [
        (base, qual)
        for base, qual in zip(record.seq, record.letter_annotations["phred_quality"])
    ]
    aligned = []
    i = 0
    for b in alignment:
        if b != "-":
            orig_base, orig_qual = original[i]
            i += 1
            aligned.append([b, orig_base, orig_qual])
        else:
            aligned.append([b, "N", missing_qual])
    seq = "".join(s for b, s, q in aligned)
    phred = [q for b, s, q in aligned]
    return seq, phred
