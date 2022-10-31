#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
from couplet.exceptions import (
    InsufficientPythonVersionError,
)

import math
import itertools
import sys
import logging
import csv
from collections import defaultdict
import pandas as pd

LOGGER = logging.getLogger("root")


def get_base_name(args):
    # generate a basname for the output files
    assert (
        "_R1" in args.fq1 and "_R2" in args.fq2
    ), "make sure the file name contains _R1 and _R2"
    base_name = args.fq1.split("_R1")[0]
    base_name2 = args.fq2.split("_R2")[0]
    assert base_name == base_name2
    return base_name


def get_outfile_names(args):

    base_name = get_base_name(args)

    out_discard_file1 = base_name + "_discarded_R1.fq.gz"
    out_discard_file2 = base_name + "_discarded_R2.fq.gz"
    resolved_NV = base_name + "_resolved.fq.gz"
    stats_file = base_name + "_couplet.yaml"
    additional_stats_file = base_name + "_couplet_additional_stats.yaml"
    return (
        out_discard_file1,
        out_discard_file2,
        resolved_NV,
        stats_file,
        additional_stats_file,
    )


def get_CEGX_colours(num_colours):

    """Function to return CEGX colour scheme.

    Given a required number of colours, this function returns a colour
    scheme that matches the CEGX colours. The CEGX colour scheme consists
    of three sets of colours: Primary, secondary and tertiary colours.
    Where possible, one should first use the primary, then the secondary,
    and finally the tertiary colours. This function obeys this principle.

    Additionally, within the above constraint, this method picks colours
    that maximises the overall difference between all the colours. This
    ensures maximum differentiability of colours.

    Args:
        num_colours: Required number of colours (maximum 13).

    Returns:
        A list of colours in hex values.

    """

    if num_colours > 13:
        raise ValueError(
            f"The requested number of colours ("
            + str(num_colours)
            + ") is too large. There are at most "
            + str(13)
            + " colours available."
        )

    # Define colours: [ hex, [R, G, B] ]
    teal = ["#63B1BC", [99, 177, 188]]
    duck_egg = ["#B5E3D8", [181, 227, 216]]
    lime = ["#97D700", [151, 215, 0]]
    denim = ["#253746", [37, 55, 70]]
    dark_blue = ["#4F758B", [79, 117, 139]]
    cement = ["#D0D0CE", [208, 208, 206]]
    red = ["#CD4747", [205, 71, 71]]
    tert1 = ["#38947E", [56, 148, 126]]
    tert2 = ["#6BC7B1", [107, 199, 177]]
    tert3 = ["#96B1C7", [150, 177, 199]]
    salmon = ["#FB9A9A", [251, 154, 154]]
    pale_orange = ["#FDBE6F", [253, 190, 111]]
    tert4 = ["#A5CCDF", [165, 204, 263]]

    # Define colour groups
    primary_colours = [teal, duck_egg, lime, denim]
    secondary_colours = [dark_blue, cement, red]
    tertiary_colours = [tert1, tert2, tert3, salmon, pale_orange, tert4]

    if num_colours == 1:
        return ["#63B1BC"]
    if num_colours < 4:
        return optimise_colours([], primary_colours, num_colours)
    if num_colours == 4:
        return [col[0] for col in primary_colours]
    if num_colours < 7:
        return optimise_colours(primary_colours, secondary_colours, num_colours)
    if num_colours == 7:
        return [col[0] for col in primary_colours + secondary_colours]
    if num_colours < 13:
        return optimise_colours(
            primary_colours + secondary_colours, tertiary_colours, num_colours
        )
    if num_colours == 13:
        return [
            col[0] for col in primary_colours + secondary_colours + tertiary_colours
        ]


def optimise_colours(fixed_colours, target_colours, num_colours):

    """Function to identify a set of colours with maximal differentiability.

    This function takes sets of fixed and target colours. The fixed set
    represents the colours which MUST be in the final output set. The fixed
    set can be empty ([]). The target set represents the colours from which
    we are to identify a subset that goes into the final output set.

    The final output set is chosen such that the identified colours have
    the maximum possible overall difference between them. For example, if

      Fixed colours: Green, blue
      Target colours: Black, purple, red, light blue

    and we are asked to pick four colours, we can pick:

       - Green, blue, black, purple
       - Green, blue, black, red
       - Green, blue, black, light blue
       - Green, blue, purple, red
       - Green, blue, purple, light blue
       - Green, blue, red, light blue

    For each of the above, the method will compute:

      D = dist(col1, col2) + dist(col1, col3) + dist(col1, col4) +
          dist(col2, col3) + dist(col2, col4) + dist(col3, col4)

    where col1, col2, col3, and col4 are the four colours. The dist() function
    is the standard Euclidean distance function operating on RGB values. The set
    of four colours with the maximal value of D is returned. In the above case,
    purple and light blue are very close to blue, so we would expect the function
    to pick black and red as the target colours.

    Colours are in the format: [ hex, [ R, G, B ] ].

    Args:
        fixed_colours: A list of colours that must be in the output set.
        target_colours: A list of colours for which a subset will be added
                        to the output set. This subset will be the set of
                        colours that maximises the overall distance between
                        the colours in the output set.
        num_colours: Required number of colours.

    Returns:
        A list of optimal colours in hex values.

    """

    if num_colours < len(fixed_colours):
        raise ValueError(
            f"The requested number of colours ("
            + str(num_colours)
            + ") is too small. At least "
            + str(len(fixed_colours))
            + " colours should be used."
        )

    if num_colours > len(fixed_colours + target_colours):
        raise ValueError(
            f"The requested number of colours ("
            + str(num_colours)
            + ") is too large. There are at most "
            + str(len(fixed_colours + target_colours))
            + " colours available."
        )

    # Define combinations of target colours
    num_pick = num_colours - len(fixed_colours)
    all_combs = itertools.combinations(target_colours, num_pick)
    total_distances = []
    colour_groups = []

    # Loop over all combinations of colours and compute the overall distance
    for comb in all_combs:

        colour_group = fixed_colours + list(comb)
        colour_groups.append(colour_group)
        dist = 0
        pair_combs = itertools.combinations(colour_group, 2)

        for pair in pair_combs:
            pair_colours = list(pair)

            # Compute Euclidean distance on RGB colours
            dist += math.dist(pair_colours[0][1], pair_colours[1][1])

        total_distances.append(dist)

    # Return the colours with the maximal overall distance between then
    best_colours = colour_groups[total_distances.index(max(total_distances))]
    return [col[0] for col in best_colours]


def check_python_version(min_version):

    """Function to ensure that a minimum Python version is used.

    Given a minimum Python version ([major, minor, micro]), the current
    Python version is checked against the minimum. If the current version
    is less than the minimum, an InsufficientPythonVersionError is thrown.

    Args:
        min_version: A minimum required Python version. Format: [ major,
        minor, micro ].

    Returns:
        [ Nothing ]
    """

    version = sys.version_info[0:3]  # Current version: [Major, minor, micro]
    if (
        version[0] < min_version[0]
        or (version[0] == min_version[0] and version[1] < min_version[1])
        or (
            version[0] == min_version[0]
            and version[1] == min_version[1]
            and version[2] < min_version[2]
        )
    ):
        raise InsufficientPythonVersionError(min_version, version)


def load_quality_table(path):
    """Function to load a phred table into a dictionary.

    Args:
        path: A path to a phred table in CSV format like:
               first_base,second_base,first_phred,second_phred, correct, incorrect, error_rate,         phred
                       T,           T,         37,          37,15296191,      2471, 0.00016151739282820943,37
               ....

    Returns:
        [ Nothing ]
    """
    quality_dict = dict()
    df = pd.read_csv(path, comment="#")
    for row in df.itertuples():
        b1 = row.first_base
        b2 = row.second_base
        q1 = row.first_phred
        q2 = row.second_phred
        quality_dict[(b1, b2, q1, q2)] = row.phred

    LOGGER.info("Quality table used:")
    for key, value in quality_dict.items():
        LOGGER.info("\t" + str(key) + " : " + str(value))

    return quality_dict
