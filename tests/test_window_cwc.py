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

    seqs = pytestconfig.rootdir / "tests" / "data" / "raw_data" / "window_cwc.fa"

    # Seqs must be converted to a str fro some reason...
    model.number(str(seqs))

    return model


def test_files_are_identical(anarcii_model, tmp_path, pytestconfig):
    expected_file_templates = {
        "txt": (
            pytestconfig.rootdir / "tests" / "data/expected_data/window_expected_1.txt"
        ),
        "json": (
            pytestconfig.rootdir / "tests" / "data/expected_data/window_expected_1.json"
        ),
    }

    # Generate and check both text and json files
    for filetype in ["txt", "json"]:
        test_file = tmp_path / f"window_test_1.{filetype}"
        expected_file = expected_file_templates[filetype]

        # if filetype == "txt":
        #     anarcii_model.to_text(test_file)
        # else:
        #     anarcii_model.to_json(test_file)

        ls_format = anarcii_model.to_list()

        # need to write to a json
        with open(test_file, "w") as f:
            json.dump(ls_format, f)

        if filetype == "json":
            with open(expected_file) as f1, open(test_file) as f2:
                json_expected = json.load(f1)
                json_test = json.load(f2)

            # Ensure both lists have the same length
            assert len(json_expected) == len(json_test), (
                f"Expected list length {len(json_expected)} but got {len(json_test)}"
            )

            # Iterate over both lists concurrently
            for expected_item, test_item in zip(json_expected, json_test):
                expected_number, expected_data = expected_item
                test_number, test_data = test_item

                assert expected_number == test_number, (
                    f"Numbering for {expected_data['query_name']} is different! "
                    f"Expected: {expected_number}, Got: {test_number}"
                )
                reference = pytest.approx(expected_data["score"], abs=0.5)
                assert test_data["score"] == reference, (
                    f"Scores differ more than 0.5 for {expected_data['query_name']}! "
                    f"Expected: {expected_data['score']}, Got: {test_data['score']}"
                )
