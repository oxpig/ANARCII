import json

import pytest

from anarcii import Anarcii


@pytest.fixture(scope="session")
def anarcii_model(pytestconfig):
    model = Anarcii(
        seq_type="antibody",
        batch_size=64,
        cpu=False,
        ncpu=12,
        mode="speed",
        verbose=False,
    )
    seqs = pytestconfig.rootpath / "tests" / "data" / "raw_data" / "sabdab_filtered.fa"

    model.number(seqs)

    return model


@pytest.mark.parametrize("scheme", ["chothia", "imgt", "martin", "kabat", "aho"])
def test_files_are_identical(anarcii_model, tmp_path, scheme, pytestconfig):
    suffix = f"_{scheme}"

    expected_file = (
        pytestconfig.rootpath
        / "tests"
        / "data"
        / "expected_data"
        / f"antibody{suffix}_expected_1.json"
    )

    test = anarcii_model.to_scheme(scheme)

    # Convert the test dict to a list of values
    test = list(test.values())

    with open(expected_file) as f1:
        expected = json.load(f1)

    # Ensure both lists have the same length
    assert len(expected) == len(test), (
        f"Expected list length {len(expected)} but got {len(expected)}"
    )

    # Iterate over both lists concurrently
    for expected_item, test_item in zip(expected, test):
        expected_number, expected_data = expected_item
        test_number, test_data = test_item["numbering"], test_item

        # The json files currently drop all tuples, so we need to undo this.
        # Tuple rebuild.
        expected_number = [((x[0][0], x[0][1]), x[1]) for x in expected_number]

        assert expected_number == test_number, (
            f"Numbering for {expected_data['query_name']} is different! "
            f"Expected: {expected_number}, Got: {test_number}"
        )
        reference = pytest.approx(expected_data["score"], abs=0.5)
        assert test_data["score"] == reference, (
            f"Scores differ more than 0.5 for {expected_data['query_name']}! "
            f"Expected: {expected_data['score']}, Got: {test_data['score']}"
        )
