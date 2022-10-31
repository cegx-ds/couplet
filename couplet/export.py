#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import yaml
import matplotlib.pyplot as plt
import seaborn as sns
import re
import pandas as pd
import sys
import logging
import collections
import math
from couplet.utils import (
    get_CEGX_colours,
    optimise_colours,
)
from couplet.exceptions import IncompatibleErrorCodes

LOGGER = logging.getLogger("root")

# Predefine frequently used variables
groups = [
    "acceptable_naive",
    "rescued_naive",
    "rescued_aligned",
    "discarded_naive",
    "discarded_aligned",
]

groups_nice_formatting = [
    "Acceptable reads: Naive alignment",
    "Rescued reads: Naive alignment",
    "Rescued reads: NW alignment",
    "Discarded reads: Naive alignment",
    "Discarded reads: NW alignment",
]


def log_core_stats(stats, stats_file):

    """Write core stats from SS-resolve to logger and yaml file.

    This function extracts number of input, resolved, perfect, alignable and alignment
    improved reads from the stats dictionary. The number of mismatches pre and post
    alignment are also logged.

    Args:
        stats: Dictionary storing statistics from SS-resolve.
        stats_file: Path to output file for storing core stats metrics.

    Returns:
        [ nothing ]
    """

    # Extract statistics from stats object and store in dictionary
    core_stats = {
        "input_reads": stats.get("num_reads", 0),
        "mismatches_naive_alignment": stats.get("num_mismatch_naive_alignment", 0),
        "acceptable_reads": stats.get("acceptable", 0),
        "potentially_rescuable_reads": stats.get("num_reads", 0)
        - stats.get("acceptable", 0),
        "mismatches_after_attempted_rescue": stats.get(
            "num_mismatch_post_attempted_rescue", 0
        ),
        "rescued_reads": stats.get("rescued", 0),
        "resolved_reads": stats.get("num_reads", 0) - stats.get("discarded", 0),
        "discarded_reads": stats.get("num_reads", 0)
        - stats.get("acceptable", 0)
        - stats.get("rescued", 0),
        "input_bases": stats.get("input_bases", 0),
        "produced_bases": stats.get("produced_bases", 0),
    }

    # Write core stats to logging object
    for k in core_stats:
        LOGGER.info("Core statistic %s: %.4f", k, core_stats[k])

    # Write core stats to file
    with open(stats_file, "w") as fp:
        yaml.dump(core_stats, fp, sort_keys=False)


def log_core_stats_merged(input_files, output_folder, output_prefix):

    """Function to summarise and log core stats from multiple data splits.

    This function takes a list of input yaml files with stats from mutliple
    data splits. The code merges this data by adding up values for dictionary
    entries with identical keys. The output is written as a yaml file with
    path {output_folder}/{output_prefix}_couplet.yaml.

    Additionally, log files for the data splits are merged into a single log
    file.

    Args:
        input_files: List of input stats files (one for each split).
        output_folder: Path to folder for storing outputs.
        output_prefix: Prefix to be appended to the output file.

    """

    stats = None
    for input_file in input_files:
        with open(input_file, "r") as fp:
            if stats is None:
                stats = yaml.safe_load(fp)
            else:
                next_stats = yaml.safe_load(fp)
                for k, v in next_stats.items():
                    stats[k] += v

    stats["average_mismatches_naive_alignment"] = (
        stats["mismatches_naive_alignment"] / stats["input_reads"]
    )

    if stats["potentially_rescuable_reads"] == 0:
        stats["average_mismatches_after_attempted_rescue"] = 0

    else:
        stats["average_mismatches_after_attempted_rescue"] = (
            stats["mismatches_after_attempted_rescue"]
            / stats["potentially_rescuable_reads"]
        )

    for k in list(stats.keys()):
        if (
            k.startswith("mismatches")
            or k.startswith("average")
            or k.startswith("input_reads")
            or k.startswith("input_bases")
        ):
            continue
        elif k == "rescued_reads":
            if stats["potentially_rescuable_reads"] == 0:
                stats[f"rescued_reads_efficiency_percentage"] = 0
            else:
                stats[f"rescued_reads_efficiency_percentage"] = (
                    stats[k] / stats["potentially_rescuable_reads"] * 100
                )
        elif k == "produced_bases":
            stats["produced_bases_percentage"] = stats[k] / stats["input_bases"] * 100
        else:
            stats[f"{k}_percentage"] = stats[k] / stats["input_reads"] * 100

    with open(output_folder + "/" + output_prefix + "_couplet.yaml", "w") as fp:
        yaml.dump(stats, fp, sort_keys=False)

    # Merge log files (file-by-file; not chonological entries)
    # Assumes input_files are of form path/SPLIT_*couplet.yaml
    # And log files are of form path/SPLIT_*couplet.log
    log_files = [file.replace(".yaml", ".log") for file in input_files]
    merged_log_filepath = output_folder + "/" + output_prefix + "_couplet.log"
    with open(merged_log_filepath, "w") as merged_log_file:
        for fname in log_files:
            with open(fname) as log_file:
                merged_log_file.write(log_file.read())


def log_additional_stats(stats, rule, stats_file):

    """Write additional stats from SS-resolve to logger and output file.

    This function extracts number and location of mismatches for the error codes
    defined in rule.error_codes. Mismatches are extracted for three sets
    of reads:

      - Acceptable reads: Reads that doesn't need pairwise alignment
      - Rescued reads: Reads that required pairwise alignment
      - Discarded reads: Reads that, despite pairwise alignment, didn't meet
        the quality requirement.

    For the latter two categories, mismatches are logged both before and after
    pairwise resolution.

    The function also logs read length information.

    Args:
        stats: Dictionary storing statistics from SS-resolve.
        rule: ResolutionRule object.
        stats_file: Path to output file for storing mismatch stats metrics.
                    Assumed to end in ".yaml".

    Returns:
        [ nothing ]
    """

    # Loop over the five groups and log statistics
    for i in range(len(groups)):

        group = groups[i]

        # Define absolute mismatch stats by extracting relevant values from the stats dictionary
        keys = [group + "_" + feature for feature in rule.error_codes]
        values = [stats.get(feature, 0) for feature in keys]
        mismatch_stats = dict(zip(keys, values))

        # Log mismatch stats (absolute and normalised) to logging object
        for j in mismatch_stats:
            LOGGER.info("Mismatch statistic %s: %.4f", j, mismatch_stats[j])

        # Write mismatch stats to output file
        if i == 0:
            with open(stats_file, "w") as fp:
                yaml.dump(mismatch_stats, fp)
        else:
            with open(stats_file, "a") as fp:
                yaml.dump(mismatch_stats, fp)

    """
    #
    # Write mismatch positional data
    #
    """
    pos_keys = [key for key in stats.keys() if re.search(rf"_pos_", key)]
    values = [stats.get(key, 0) for key in pos_keys]
    mismatch_pos_stats = dict(zip(pos_keys, values))
    if len(mismatch_pos_stats) > 0:
        with open(stats_file, "a") as fp:
            yaml.dump(mismatch_pos_stats, fp)

    # Log the error code for use downstream
    with open(stats_file, "a") as fp:
        yaml.dump({"error_codes_used_" + ",".join(rule.error_codes): 0}, fp)

    """
    #
    # Write read length data
    #
    """
    length_keys = [key for key in stats.keys() if re.search(rf"_length_", key)]
    length_values = [stats.get(key, 0) for key in length_keys]
    length_stats = dict(zip(length_keys, length_values))
    with open(stats_file, "a") as fp:
        yaml.dump(length_stats, fp)


def log_and_plot_additional_stats(stats, group_numbers, error_codes, stats_file):

    """Log and plot additional stats from SS-resolve.

    Given a stats dictionary and a set of error codes, this function extracts
    mismatch statistics (number and location of mismatches) and read length
    distribution metrics. The stats are logged and plots are produced.
    Mismatches are extracted for three sets of reads:

      - Acceptable reads: Reads that doesn't need pairwise alignment
      - Rescued reads: Reads that required pairwise alignment
      - Discarded reads: Reads that, despite pairwise alignment, didn't meet
        the quality requirement.

    For the latter two categories, mismatches are logged both before and after
    pairwise resolution.

    Two plots are generated. First, a bar plot is produced showing the number
    of different  mismatch types for each group of reads. Second, a lines plot
    is produced showing the distribution of these mismatches across read cycles.


    Args:
        stats: Dictionary storing statistics from SS-resolve.
        group_numbers: Number of reads in each of the read groups.
        error_codes: A list of error codes, e.g. ["0","1",...,"n"].
        stats_file: Path to output file for storing additional stats metrics.
                    Assumed to end in "_additional_stats.yaml".

    Returns:
        [ nothing ]
    """

    # Loop over the five groups and log/plot statistics
    for i in range(len(groups)):

        group = groups[i]
        group_number = group_numbers[i]
        if group_number == 0:
            continue

        # Define absolute mismatch stats by extracting relevant values from the stats dictionary
        keys = [group + "_" + feature for feature in error_codes]
        values = [stats.get(feature, 0) for feature in keys]
        mismatch_stats = dict(zip(keys, values))

        # Define normalised mismatch stats by dividing by the number of reads in the current group
        keys_norm = [key + "_norm" for key in keys]
        values_norm = [val / group_number for val in values]
        mismatch_stats_norm = dict(zip(keys_norm, values_norm))

        # Log mismatch stats (absolute and normalised) to logging object
        for j in mismatch_stats:
            LOGGER.info("Mismatch statistic %s: %.4f", j, mismatch_stats[j])
        for j in mismatch_stats_norm:
            LOGGER.info("Mismatch statistic %s: %.4f", j, mismatch_stats_norm[j])

        # Write mismatch stats to output file
        if i == 0:
            with open(stats_file, "w") as fp:
                yaml.dump(mismatch_stats, fp, sort_keys=False)
                yaml.dump(mismatch_stats_norm, fp, sort_keys=False)
        else:
            with open(stats_file, "a") as fp:
                yaml.dump(mismatch_stats, fp, sort_keys=False)
                yaml.dump(mismatch_stats_norm, fp, sort_keys=False)

        """
        #
        # Generate (normalised) mismatch counts plots
        #
        """
        # Create barplot with bars = normalised number of mismatches for each mismatch code
        df = pd.DataFrame()
        df["Mismatch code"] = error_codes
        df["Average per-read occurrence"] = values_norm
        ax = sns.barplot(
            x="Mismatch code", y="Average per-read occurrence", data=df, color="#63B1BC"
        )
        ax.set_title(groups_nice_formatting[i])

        # Define output filename for plot and save as png file
        plot_name = stats_file.replace(
            "_additional_stats.yaml", "_" + group + "_mismatch_stats_counts.png"
        )
        plt.savefig(plot_name)
        plt.clf()

        """
        #
        # Generate (normalised) mismatch position plots
        #
        """
        # Identify the maximum position logged across all the error codes
        # We do this by searching for stats keys which match "_pos_". We then extract the pos from the key.
        max_pos = int(
            max(
                [
                    int(key.split("_pos_", 1)[1])
                    for key in stats.keys()
                    if re.search(rf"_pos_", key)
                ]
            )
        )

        # Create dataframe to hold data
        df = pd.DataFrame()
        df["Read cycle"] = range(max_pos + 1)

        # Loop over mismatch codes
        for code in error_codes:

            # Define key names for postions 0,1,2,...,max_pos
            # Key names are e.g. "pwa_aligned_3_pos_88" where "88" is the position
            pos_keys = [
                group + "_" + code + "_pos_" + str(pos) for pos in range(max_pos + 1)
            ]

            # Extract values for the key names. Use a value of zero if key name is not found in stats.
            # Normalise the values by group size.
            values = [stats.get(key, 0) / group_number for key in pos_keys]

            # Add values to data frame
            df[code] = values

        # Format plot
        df = df.melt(
            "Read cycle",
            var_name="Mismatch code",
            value_name="Average per-read occurrence",
        )
        cegx_colours = get_CEGX_colours(len(error_codes))
        ax = sns.lineplot(
            data=df,
            x="Read cycle",
            y="Average per-read occurrence",
            hue="Mismatch code",
            palette=cegx_colours,
        )
        ax.set_title(groups_nice_formatting[i])

        # Shrink current axis by 14%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.86, box.height])

        # Put a legend to the right of the current axis
        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), title="Mismatch code")

        # Define output filename for plot and save as png file
        plot_name = stats_file.replace(
            "_additional_stats.yaml", "_" + group + "_mismatch_stats_pos.png"
        )
        plt.savefig(plot_name)
        plt.clf()

    # Log read length distribution metrics
    # Get min, mean, median and max lengths
    (
        acceptable_L_min,
        acceptable_L_mean,
        acceptable_L_std_dev,
        acceptable_L_median,
        acceptable_L_max,
    ) = get_length_metrics("acceptable", stats)
    (
        rescued_L_min,
        rescued_L_mean,
        rescued_L_std_dev,
        rescued_L_median,
        rescued_L_max,
    ) = get_length_metrics("rescued", stats)
    (
        discarded_L_min,
        discarded_L_mean,
        discarded_L_std_dev,
        discarded_L_median,
        discarded_L_max,
    ) = get_length_metrics("discarded", stats)

    # Add to length stats
    length_stats = {"acceptable_min_length": acceptable_L_min}
    length_stats["acceptable_mean_length"] = acceptable_L_mean
    length_stats["acceptable_std_dev_length"] = acceptable_L_std_dev
    length_stats["acceptable_median_length"] = acceptable_L_median
    length_stats["acceptable_max_length"] = acceptable_L_max

    length_stats["rescued_min_length"] = rescued_L_min
    length_stats["rescued_mean_length"] = rescued_L_mean
    length_stats["rescued_std_dev_length"] = rescued_L_std_dev
    length_stats["rescued_median_length"] = rescued_L_median
    length_stats["rescued_max_length"] = rescued_L_max

    length_stats["discarded_min_length"] = discarded_L_min
    length_stats["discarded_mean_length"] = discarded_L_mean
    length_stats["discarded_std_dev_length"] = discarded_L_std_dev
    length_stats["discarded_median_length"] = discarded_L_median
    length_stats["discarded_max_length"] = discarded_L_max

    # Log mismatch stats (absolute and normalised) to logging object
    for j in length_stats:
        LOGGER.info("Length statistic %s: %.4f", j, length_stats[j])

    with open(stats_file, "a") as fp:
        yaml.dump(length_stats, fp, sort_keys=False)


def log_and_plot_additional_stats_single(stats, rule, stats_file):

    """Orchestration function for logging additional stats for a single
    pair of reads files.

    This function is intended for use in cases where a single pair of
    reads files (name_R1.fq.gz, name_R2.fq.gz) are analysed. This is
    in contrast to parallel analysis where the pair of read files are
    split into smaller chunks.

    Mismatches are extracted for three sets of reads:

      - Acceptable reads: Reads that doesn't need pairwise alignment
      - Rescued reads: Reads that required pairwise alignment
      - Discarded reads: Reads that, despite pairwise alignment, didn't meet
        the quality requirement.

    Read length distribution metrics are computed.

    This function defines the number of reads in each category, then calls
    a secondary function responsible for logging and plotting the mismatch
    data.

    Args:
        stats: Dictionary storing statistics from SS-resolve.
        rule: ResolutionRule object.
        stats_file: Path to output file for storing additional stats metrics.
                    Assumed to end in ".yaml".

    Returns:
        [ nothing ]
    """

    # Get number of reads in the different groups of reads
    num_acceptable = stats.get("acceptable", 0)
    num_rescued = stats.get("rescued", 0)
    num_discarded = stats.get("discarded", 0)

    # Define number of reads in each group
    group_numbers = [
        num_acceptable,
        num_rescued,
        num_rescued,
        num_discarded,
        num_discarded,
    ]

    log_and_plot_additional_stats(stats, group_numbers, rule.error_codes, stats_file)


def log_and_plot_additional_stats_merged(input_files, output_folder, output_prefix):

    """Function to summarise, log and plot additional stats from multiple data splits.

    This function is intended for use in cases where, due to a parallel
    analysis approach, the read files have been split into smaller chunks.

    Mismatches are extracted for three sets of reads:

      - Acceptable reads: Reads that doesn't need pairwise alignment
      - Rescued reads: Reads that required pairwise alignment
      - Discarded reads: Reads that, despite pairwise alignment, didn't meet
        the quality requirement.

    Read length distribution metrics are logged.

    This function defines the number of reads in each category, then calls
    a secondary function responsible for logging and plotting the mismatch
    data.

    This function assumes that 'log_core_stats_merged' has already been run.

    Args:
        input_files: List of input mismatch stats files (one for each split).
        output_folder: Path to folder for storing outputs.
        output_prefix: Prefix to be appended to the output files.
    """

    # Create a "merged" stats dictionary by adding entries from all splits
    stats = None
    for input_file in input_files:
        with open(input_file, "r") as fp:
            if stats is None:
                stats = yaml.safe_load(fp)
            else:
                next_stats = yaml.safe_load(fp)
                for k, v in next_stats.items():
                    if k in stats:
                        stats[k] += v
                    else:
                        stats[k] = v

    # Load total number of reads in each category from core stats file
    with open(output_folder + "/" + output_prefix + "_couplet.yaml", "r") as fp:
        core_stats = yaml.safe_load(fp)

    # Get number of reads in the different groups of reads
    num_acceptable = core_stats.get("acceptable_reads", 0)
    num_rescued = core_stats.get("rescued_reads", 0)
    num_discarded = core_stats.get("discarded_reads", 0)

    # Define number of reads in each group
    group_numbers = [
        num_acceptable,
        num_rescued,
        num_rescued,
        num_discarded,
        num_discarded,
    ]

    # Extract the error code from a special entryin the stats dict.
    # The entry is of form "rescued_aligned_error_codes_used_0,1,2,...,n".
    # In the following we extract just the error codes ("0,1,2,...,n").
    # Upon creating a set, there should only be one set of error codes left.
    # If not, this tells us that the input (split) files were incompatbile.
    error_codes_string = set(
        [
            key.split("_used_")[1]
            for key in stats.keys()
            if re.search(rf"error_codes_used_", key)
        ]
    )
    if len(error_codes_string) > 1:

        raise IncompatibleErrorCodes(error_codes_string)

    error_codes = list(error_codes_string)[0].split(",")

    additional_stats_file = (
        output_folder + "/" + output_prefix + "_couplet_additional_stats.yaml"
    )
    log_and_plot_additional_stats(
        stats, group_numbers, error_codes, additional_stats_file
    )


def get_length_metrics(name, stats):

    """Function to compute stats on read length distribution.

    Given a "stats" dictionary, this function computes the minimum,
    mean, (population) standard deviation, median and maximum read
    length. A name input is used to define the type of reads to
    compute stats on, e.g. "acceptable", "rescued", "discarded".

    The function purposefully works on counts data, e.g.
    {"acceptable_length_85: 1021", "acceptable_length_86": 1283, ...}.
    This is to avoid having to convert the counts into overly long
    lists, e.g. [85, 85, ..., 85, 86, 86, ..., 86, ..., 122].

    Note that the standard deviation is the population standard
    deviation, not the sample standard deviation.

    Args:
        name: A qualifier used in the stats dictionary keys. The
            dictionary keys are assumed to be of the form:
            "{name}_length_{read_length}".
        stats: A dictionary object containing information about
            read lengths.

    Returns
        min, mean, std_dev, median and max read lengths

    """

    # Create counts dictionary of form {1:5, 2:10, 3:8, 4:12, ..}
    counts_dict = {
        int(key.split(f"{name}_length_")[1]): value
        for key, value in stats.items()
        if re.search(rf"{name}_length_", key)
    }

    # Sort by key (i.e. length of read)
    sorted_dict = collections.OrderedDict(sorted(counts_dict.items()))

    # Compute mean: (1*5 + 2*10 + 3*8 + 4*12 +...) / (5+10+8+12+...)
    weighted_counts = sum(length * count for length, count in sorted_dict.items())
    total_counts = sum(sorted_dict.values())
    mean_L = weighted_counts / total_counts

    # Compute standard deviation
    total_squares = sum(
        length * length * count for length, count in sorted_dict.items()
    )
    mean_of_squares = total_squares / total_counts
    variance = mean_of_squares - mean_L * mean_L
    std_dev_L = math.sqrt(variance)

    # Compute min, max
    min_L = list(sorted_dict.keys())[0]
    max_L = list(sorted_dict.keys())[-1]

    # Compute median
    if total_counts % 2 == 0:
        even = True
        half_counts = total_counts / 2
    else:
        even = False
        half_counts = math.ceil(total_counts / 2)

    count_sum = 0
    median_L = -1
    for i in range(len(sorted_dict)):

        length = list(sorted_dict.keys())[i]
        count = list(sorted_dict.values())[i]
        count_sum += count

        if count_sum > half_counts:

            median_L = length
            break

        if count_sum == half_counts:

            if even == True:
                next_length = list(sorted_dict.keys())[i + 1]
                median_L = (length + next_length) / 2
            else:
                median_L = length
            break

    return min_L, mean_L, std_dev_L, median_L, max_L
