"""Tests for :func:`bcmproteomics.data_loading.fmt_input_data`."""

from pathlib import Path
import sys

import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bcmproteomics.data_loading import fmt_input_data


@pytest.mark.parametrize(
    "input_data, expected",
    [
        pytest.param([12345], [{"recno": 12345}], id="integer"),
        pytest.param(
            [(12346, 2)],
            [{"recno": 12346, "runno": 2, "searchno": 1}],
            id="rec-run tuple",
        ),
        pytest.param(
            [(12347, 3, 4)],
            [{"recno": 12347, "runno": 3, "searchno": 4}],
            id="rec-run-search tuple",
        ),
        pytest.param(
            [{"recno": 12348, "runno": 5, "searchno": 6}],
            [{"recno": 12348, "runno": 5, "searchno": 6}],
            id="dictionary",
        ),
        pytest.param(
            [(12349, 1, 2, 3)],
            [],
            id="invalid tuple",
        ),
    ],
)
def test_fmt_input_data_parametrized_cases(input_data, expected):
    """Fmt_input_data should standardize dictionaries and skip invalid inputs."""
    assert fmt_input_data(input_data) == expected


def test_fmt_input_data_mixed_inputs_skips_only_invalid():
    """Ensure only the invalid entries are removed from the result."""
    input_data = [
        (11111, 2),
        (22222, 3, 4, 5),
        33333,
    ]
    expected = [
        {"recno": 11111, "runno": 2, "searchno": 1},
        {"recno": 33333},
    ]
    assert fmt_input_data(input_data) == expected
