#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import pytest
from pytest import approx
import random
import collections
import statistics

from ssresolve.export import (
    get_length_metrics,
)


@pytest.mark.parametrize(
    "lengths",
    [
        [70, 78, 82, 84, 85, 85, 86, 86, 86, 86, 87, 89, 90, 94, 96],
        [5, 8, 9, 12],  # Even
        [8, 99, 120, 154],  # Odd
        [8, 8, 8, 8],  # All the same, even
        [8, 8, 8, 8, 8],  # All the same, odd
        [56],  # Single value
    ],
)
def test_get_length_metrics(lengths):

    # Shuffle lengths (with seed) so as to ensure order of
    # lengths doesn't matter
    random.Random(1).shuffle(lengths)

    # Create a counts dictionary
    counts = collections.Counter(lengths)

    # Prefix keys with some name, e.g. "rescued"
    counts = {f"rescued_length_{k}": v for k, v in counts.items()}

    print(counts)

    # Compute results
    min_L, mean_L, std_dev_L, median_L, max_L = get_length_metrics("rescued", counts)

    assert min(lengths) == min_L
    assert statistics.mean(lengths) == mean_L
    assert statistics.pstdev(lengths) == approx(std_dev_L)
    assert statistics.median(lengths) == median_L
    assert max(lengths) == max_L
