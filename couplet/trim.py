#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from couplet.align import count_mismatches, get_aligned_record


def get_dynamic_right_trimming_point(
    mismatches, read_length, min_read_length, left_trim, threshold=0.05
):

    """Function to dynamically identify a right trimming point for a read.

    Given a minimum read length and a mismatch density threshold, this function
    identifies the right trimming point that rescues the read, i.e. the point where,
    if the read is trimmed from the right, the resulting read has a mismatch
    density below the threshold. Mismatch densities are computed as "number  of
    mismatches" / "read length". The function assumes that right trailing N's have
    already been removed from the read / the list of mismatches. Additionally, the
    function assumes that a left trailing N trimming point has already been
    identified (but not applied). The code then ignores N's accounted for by the
    left trimming.


    If no adequate trimming point is identified, the code returns None.

    Args:
        mismatches: A list of mismatch indices, i.e. the location of mismatches
            in the read. This list is zero based. Right trailing mismatches
            should have already been removed.
        read_length: The length of the read after right trimming.
        min_read_length: The minimum required read length.
        left_trim: The trimming point associated with left trailing N's. E.g.
            if the sequence is "NNNACGT", then left_trim = 3. left_trim can
            take a value of zero (if there are no left trailing N's).
        threshold: The maximum mismatch density a read must have.

    Returns:
        The optimal trimming point to rescue the read. If no solution is found,
        the code returns the input read_length as the trimming point.

    """

    if read_length == 0:
        return 0

    # Trim left trailing N's from list of mismatches, e.g.
    # [0,1,2,3,10,20,30,40] --> [10,20,30,40]
    mismatches = mismatches[left_trim:]

    # Check if mismatch trimming is necessary
    if len(mismatches) <= (read_length - left_trim) * threshold:
        return read_length

    for i in range(len(mismatches) - 1, -1, -1):

        # Length after trimming
        length = mismatches[i] - left_trim

        # If length is less than requirement return
        # min_read_length, i.e. no trimming
        if length < min_read_length:
            return read_length

        # Check if threshold is met
        if i <= length * threshold:

            return mismatches[i]

    # No solution found. This should never happen.
    assert (
        False
    ), "Unreachable statement in couplet.trim.get_dynamic_right_trimming_point."


def get_left_trimming_point(mismatches):

    """Function to identify a left trimming point for removing
    trailing mismatches.

    Given a list of mismatches, this function identifies a trimming
    point, where if the read was trimmed from the left, this would
    result in the removal of trailing mismatches.

    Args:
        mismatches: A list of mismatch indices, i.e. the location of mismatches
            in the read. This list is zero based

    Returns:
        The left trimming point. For example, a read "NNNACGTNAC" with
        mismatches=[0,1,2,7] would result in a left trimming point of 3.
        The trimming point can be zero.

    """

    if len(mismatches) == 0:
        return 0

    # Find trimming point from left.
    # E.g. if mismatches are [0,1,2,3,10,20] then left_trim = 4.
    left_trim = mismatches[-1] + 1
    for i, pos in enumerate(mismatches):
        if pos > i:
            left_trim = i
            break

    return left_trim


def get_right_trimming_point(mismatches, read_length):

    """Function to identify a right trimming point for removing
    trailing mismatches.

    Given a list of mismatches and a read length, this function identifies a
    trimming point, where if the read was trimmed from the right, this would
    result in the removal of trailing mismatches.

    Args:
        mismatches: A list of mismatch indices, i.e. the location of mismatches
            in the read. This list is zero based
        read_length: The length of the read to be trimmed.

    Returns:
        The right trimming point. For example, a read "ACGTNACGTNN" with
        read_length=11 and mismatches=[4,9,10] would result in a right
        trimming point of 9. The trimming point can be zero (but only in
        exceptional circumstances where all bases in the read are mismatches).

    """

    # Start with right_trim = read_length, then reverse through mismatches and
    # decrease right_trim as long as there is a succession of mismatches.
    right_trim = read_length
    for i in range(len(mismatches)):
        pos = mismatches[-(i + 1)]
        if pos < (read_length - i - 1):
            break
        else:
            right_trim = pos

    return right_trim


def trim_read(read_id, seq, phred, left_trim, right_trim):

    """Function to trim a read given left and right trimming points.

    Given trimming boundaries, this function trims the sequence and
    Phred scores of a read from the left and the right. The function
    returns the trimmed data as a SeqRecord object.

    Args:
        read_id: A read name to be added to the resulting SeqRecord.
        seq: The sequence to be trimmed. This is a string, e.g. "ACGT".
        phred: The Phred scores to be trimmed. This is a list of
            numerical values, e.g. [ 30, 31, 29, 32 ].
        left_trim: The trimming point when trimming from the left.
            Objects are trimmed using slicing FROM the left_trim value.
            As an example, left_trim=3 means trimming a read
            "NNNACGTNN" at the vertical bar: "---|ACGTNN".
        right_trim: The trimming point when trimming from the right.
            Objects are trimmed using slicing TO the right_trim value.
            As an example, right_trim=7 means trimming a read
            "NNNACGTNN" at the vertical bar: "NNNACGT|NN".

    Returns:
        A basic SeqRecord object with id, seq and Phred scores.
    """

    # Trim sequence amd phred
    seq = seq[left_trim:right_trim]
    phred = phred[left_trim:right_trim]

    # Create SeqRecord
    rec = SeqRecord(
        seq=Seq(seq),
        id=f"{read_id}",
        letter_annotations={"phred_quality": phred},
        description="",
    )

    return rec
