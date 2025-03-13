import pathlib

import pytest

from anarcii.input_data_processing import coerce_input, file_input

cases = {
    "single-sequence-string": (
        "test",
        {"sequence": "test"},
    ),
    "single-name-sequence-tuple": (
        ("name", "test"),
        {"name": "test"},
    ),
    "list-single-sequence-string": (
        ["test"],
        {"sequence-1": "test"},
    ),
    "list-multiple-sequence-strings": (
        ["test-1", "test-2"],
        {"sequence-1": "test-1", "sequence-2": "test-2"},
    ),
    "list-single-name-sequence-tuple": (
        [("name", "test")],
        {"name": "test"},
    ),
    "list-multiple-name-sequence-tuples": (
        [("name-1", "test-1"), ("name-2", "test-2")],
        {"name-1": "test-1", "name-2": "test-2"},
    ),
    "dict-single-sequence": (
        {"name": "test"},
        {"name": "test"},
    ),
    "dict-multiple-sequences": (
        {"name-1": "test-1", "name-2": "test-2"},
        {"name-1": "test-1", "name-2": "test-2"},
    ),
}


@pytest.mark.parametrize("input_data, expected", cases.values(), ids=cases.keys())
def test_coerce_input(input_data, expected):
    assert coerce_input(input_data) == expected


@pytest.mark.xfail
def test_file_input():
    assert file_input(pathlib.Path())
