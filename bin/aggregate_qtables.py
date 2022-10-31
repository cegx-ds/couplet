#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import argparse
import logging
import os

import pandas as pd

from couplet.qtables import aggregate_qtables


def main():
    parser = argparse.ArgumentParser(
        "Build a consensus quality table by averaging frequencies across multiple quality tables"
    )
    parser.add_argument("--qtables", nargs="+", required=True)
    parser.add_argument(
        "--output-dir",
        help="directory in which to write any output files",
        required=True,
    )
    parser.add_argument(
        "--tag", default="", help="A tag to use to build the output file name"
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    log_file = os.path.join(args.output_dir, f"{args.tag}_aggregate_qtable.log")
    logging.basicConfig(
        filename=log_file,
        format="%(asctime)-15s %(filename)s:%(lineno)d %(message)s",
        level=logging.INFO,
    )

    qtables = []
    header_lines = ["# Aggregated Q-table"]
    for i, filename in enumerate(args.qtables):
        qtables.append(pd.read_csv(filename, comment="#"))
        with open(filename, "r") as in_file:
            header_lines.append(f"# Q-table {i}: {os.path.basename(filename)}")
            header_lines.extend([_l.strip() for _l in in_file if _l.startswith("#")])

    result = aggregate_qtables(qtables)

    output_path = os.path.join(args.output_dir, f"{args.tag}_combined_qtable.csv")
    logging.info(f"Saving to %s", output_path)

    with open(output_path, "w") as out_file:
        out_file.write("\n".join(header_lines) + "\n")
    result.to_csv(output_path, index=False, mode="a")


if __name__ == "__main__":
    main()
