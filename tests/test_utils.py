#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
import pytest

from ssresolve.utils import (
    get_CEGX_colours,
    optimise_colours,
    check_python_version,
    load_quality_table,
)
from ssresolve.exceptions import (
    InsufficientPythonVersionError,
)


def test_get_CEGX_colours():

    # One colour should return teal
    # Four colours should return all primary colours
    # Six colours should choose cement and red over dark blue
    assert get_CEGX_colours(1) == ["#63B1BC"]
    assert get_CEGX_colours(4) == ["#63B1BC", "#B5E3D8", "#97D700", "#253746"]
    assert get_CEGX_colours(6) == [
        "#63B1BC",
        "#B5E3D8",
        "#97D700",
        "#253746",
        "#D0D0CE",
        "#CD4747",
    ]

    # Check exception is raised if more than 13 colours are requested
    with pytest.raises(ValueError):
        get_CEGX_colours(14)


def test_optimise_colours():

    # Pick four colours from six. Two colours (green and blue) must be in the output set.
    green = ["#16c410", [22, 196, 16]]
    blue = ["#236ee8", [35, 110, 232]]
    black = ["#000000", [0, 0, 0]]
    purple = ["#a11dc2", [161, 29, 194]]
    red = ["#bd092d", [189, 9, 45]]
    light_blue = ["#2ab7eb", [42, 183, 235]]

    fixed = [green, blue]
    target = [black, purple, red, light_blue]

    # Expectation: Purple and lightblue are similar to blue, which is in the fixed set.
    # Thus we expect the function to pick black and red from the set of target colours.
    expected = [green, blue, black, red]

    assert optimise_colours(fixed, target, 4) == [col[0] for col in expected]

    # Check exception is raised when requesting too few/many colours
    with pytest.raises(ValueError):
        optimise_colours(fixed, target, 1)
    with pytest.raises(ValueError):
        optimise_colours(fixed, target, 7)


def test_check_python_version():

    with pytest.raises(InsufficientPythonVersionError):
        check_python_version([99, 99, 99])

    with pytest.raises(InsufficientPythonVersionError):
        check_python_version([99, 0, 99])

    with pytest.raises(InsufficientPythonVersionError):
        check_python_version([99, 99, 0])


@pytest.fixture
def qtable_path():
    # Table contents:
    # first_base,second_base,first_phred,second_phred,correct,incorrect,error_rate,phred
    # T,T,37,37,15296191,2471,0.00016151739282820943,37
    # T,C,37,37,11498964,5110,0.0004441904667859404,33
    return "tests/qtable_test.csv"


def test_load_quality_table(qtable_path):
    # print("00dsldk")
    t = load_quality_table(qtable_path)
    d = dict()
    d[("T", "T", 37, 37)] = 37
    d[("T", "C", 37, 37)] = 33
    assert t == d
