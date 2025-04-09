import json
import subprocess

import pytest


@pytest.fixture
def run_cli():
    # ruff: noqa: S603
    def _run(args):
        return subprocess.run(
            args,
            capture_output=True,
            text=True,
        )

    return _run


@pytest.mark.parametrize(
    "args, expected",
    [
        (["anarcii", "--help"], "Specify"),
        (
            [
                "anarcii",
                "GGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGG"
                "STYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRSTDSSSWYYYYYG",
            ],
            "Error: None",
        ),
        (
            [
                "anarcii",
                "-t",
                "tcr",
                "GGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGG"
                "STYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRSTDSSSWYYYYYG",
            ],
            "Score: 28",
        ),
        (
            [
                "anarcii",
                "-t",
                "unknown",
                "GGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGG"
                "STYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRSTDSSSWYYYYYG",
            ],
            "Chain: H",
        ),
    ],
)
def test_cli_commands(run_cli, args, expected):
    """
    Basic test of the command line interface (CLI) of the Anarcii tool.
    This test checks that an output comes from the CLI for various input arguments.
    It then searches the stdout for some basic expected output, like no error, score or
    chain type etc.
    The test cases include:
    - Displaying the help message
    - Providing a sequence without specifying the type, expect no error
    - Providing a TCR sequence and checking the score starts with 28
    - Providing an unknown type and checking the chain is H
    """
    result = run_cli(args)
    assert result.returncode == 0
    assert expected in result.stdout


@pytest.mark.parametrize(
    "args",
    (
        ["anarcii", "-t", "antibody", "-m", "speed", "-o", "cli_a.json", "unknown.fa"],
        ["anarcii", "-t", "tcr", "-m", "speed", "-o", "cli_b.json", "unknown.fa"],
        ["anarcii", "-t", "shark", "-m", "speed", "-o", "cli_c.json", "unknown.fa"],
        ["anarcii", "-t", "unknown", "-m", "speed", "-o", "cli_d.json", "unknown.fa"],
    ),
    ids=("antibody", "tcr", "shark", "unknown"),
)
def test_cli_files(run_cli, args, pytestconfig, tmp_path):
    """
    More detailed tests of the CLI using the unknown fasta file.

    This test checks that the output files are generated correctly and checks against
    specific expected values in the json files.

    """
    # Copy the args list so we can modify it
    new_args = args.copy()

    # Replace the output file with one in the tmp directory
    if "-o" in new_args:
        idx = new_args.index("-o")
        output_file = tmp_path / new_args[idx + 1]
        new_args[idx + 1] = str(output_file)

    # Replace the raw input file with the proper file in the raw_data folder
    new_args[-1] = str(
        pytestconfig.rootdir / "tests" / "data" / "raw_data" / "unknown.fa"
    )

    result = run_cli(new_args)
    assert result.returncode == 0

    # Load the generated JSON file from the tmp directory
    with open(output_file) as f_test:
        json_test = json.load(f_test)

    expected_file = (
        pytestconfig.rootdir / "tests" / "data" / "expected_data" / output_file.name
    )

    with open(expected_file) as f_exp:
        json_expected = json.load(f_exp)

    assert len(json_test) == len(json_expected), (
        f"Expected list length {len(json_expected)} but got {len(json_test)}"
    )

    # Compare each item (assuming the JSON is a list of [number, data] pairs)
    for expected_item, test_item in zip(json_expected, json_test):
        expected_number, expected_data = expected_item
        test_number, test_data = test_item

        query_name = expected_data["query_name"]
        assert expected_number == test_number, (
            f"Numbering for {query_name} is different! "
            f"Expected: {expected_number}, Got: {test_number}"
        )
        reference = pytest.approx(expected_data["score"], abs=0.5)
        assert test_data["score"] == reference, (
            f"Scores differ more than 0.5 for {query_name}! "
            f"Expected: {expected_data['score']}, Got: {test_data['score']}"
        )
