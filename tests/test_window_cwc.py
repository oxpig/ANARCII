import json

import pytest

from anarcii import Anarcii


@pytest.fixture(scope="session")
def anarcii_model(pytestconfig):
    model = Anarcii(
        seq_type="antibody",
        batch_size=1,
        cpu=False,
        ncpu=4,
        mode="accuracy",
        verbose=True,
    )

    seqs = pytestconfig.rootpath / "tests" / "data" / "raw_data" / "window_cwc.fa"

    results = model.number(seqs)

    return model, results


def test_files_are_identical(anarcii_model, tmp_path, pytestconfig):
    expected_file = (
        pytestconfig.rootpath / "tests" / "data/expected_data/window_expected_1.json"
    )

    _, results = anarcii_model  # Unpack results from the fixture

    # Convert the results dict to a list of values
    test = list(results.values())

    with open(expected_file) as f1:
        expected = json.load(f1)

    # Ensure both lists have the same length
    assert len(expected) == len(test), (
        f"Expected list length {len(expected)} but got {len(test)}"
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
