#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import itertools as it
import numpy as np

CHARLIST = ["A", "C", "G", "T", "N", "-"]

five_bp_rule = {
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
    # TODO create rules for each gap
}

five_bp_error_codes = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "n"]
five_bp_modifications = {"P": ["C", "mh"]}

six_bp_rule = {
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
    ("G", "G"): ("G", "Q"),
    ("G", "T"): ("N", "7"),
    ("T", "A"): ("N", "8"),
    ("T", "C"): ("C", "C"),
    ("T", "G"): ("N", "9"),
    ("T", "T"): ("T", "T"),
}

six_bp_error_codes = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
six_bp_modifications = {"PQ": ["C", "m"], "PG": ["C", "h"], "P": ["C", "mh"]}

# Q: G protected on copy strand
# PQ: mC
# CG: C
# PG: 5hmC
# CQ: Unexpected


def get_permitted_pairs(base_rule):
    """
    Returns a set of keys from a base_rule when there is no N in the genomic context values
    """
    return {k for k, v in base_rule.items() if "N" not in v[0]}


# A function for determining the score between any two bases in alignment
def create_nw_score_table(
    permitted_pairs, char_list, match_award, gap_penalty, mismatch_penalty
):
    """Create the Needleman-Wunsch score table for a given base resolution rule

    This builds a lookup table providing the Needleman-Wunsch score
    for each pair of characters from the alphabet ['A', 'C', 'G', 'T',
    'N', '-'].

    Args:
        permitted_pairs (iterable): the pairs of characters which are permitted by
            the rule
        match_award (int): the penalty applied when the pair is a match or when one
            of the two characters is N
        gap_penalty (int): the penalty applied when either of the characters in the
            pair is a '-'
        mismatch_penalty (int): the penalty applied when the pair correspond to a pair of
            valid characters (i.e., not '-' or 'N') which are not permitted under the rule

    Returns:
        the score table
    """
    k = len(char_list)
    score_table = np.zeros((k, k))
    for i, j in it.product(range(k), range(k)):
        a, b = char_list[i], char_list[j]
        if a == "-" or b == "-":
            score_table[i, j] = gap_penalty
        elif (a, b) in permitted_pairs or "N" in [a, b]:
            score_table[i, j] = match_award
        else:
            score_table[i, j] = mismatch_penalty
    return score_table


class ResolutionRule:
    def __init__(
        self,
        rule_dict,
        error_codes,
        modifications,
        match_award,
        gap_penalty,
        mismatch_penalty,
        end_gap_penalty,
    ):
        self.rule_dict = rule_dict
        self.error_codes = error_codes
        self.modifications = modifications
        self.permitted_pairs = get_permitted_pairs(rule_dict)
        self.char_list = CHARLIST
        self.index = dict(zip(self.char_list, range(len(self.char_list))))
        self.nw_score_table = create_nw_score_table(
            self.permitted_pairs,
            self.char_list,
            match_award=match_award,
            gap_penalty=gap_penalty,
            mismatch_penalty=mismatch_penalty,
        )
        self.match_award = match_award
        self.gap_penalty = gap_penalty
        self.mismatch_penalty = mismatch_penalty
        self.end_gap_penalty = end_gap_penalty
