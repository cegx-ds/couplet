#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import argparse
import yaml
import logging
from couplet.export import (
    log_core_stats_merged,
    log_and_plot_additional_stats_merged,
)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input-stats-files", nargs="+", required=True)
    parser.add_argument("--input-additional-stats-files", nargs="+")
    parser.add_argument("--output-folder", default=".")
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    # Setup LOGGER
    log_output_file = (
        args.output_folder + "/" + args.output_prefix + "_post_process.log"
    )
    logging.basicConfig(
        filename=log_output_file,
        format="%(asctime)-15s %(filename)s:%(lineno)d %(message)s",
    )
    LOGGER = logging.getLogger("root")
    LOGGER.setLevel(logging.INFO)

    # Standard stats
    log_core_stats_merged(
        args.input_stats_files, args.output_folder, args.output_prefix
    )

    # Mismatch stats
    if args.input_additional_stats_files is not None:
        log_and_plot_additional_stats_merged(
            args.input_additional_stats_files, args.output_folder, args.output_prefix
        )


if __name__ == "__main__":
    main()
