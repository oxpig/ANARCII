import pytest

from anarcii import Anarcii
from anarcii.utils import from_msgpack_map


@pytest.fixture(scope="session")
def anarcii_model(pytestconfig):
    model = Anarcii(
        seq_type="unknown",
        batch_size=16,
        cpu=False,
        ncpu=12,
        mode="speed",
        verbose=False,
    )
    seqs = pytestconfig.rootpath / "tests" / "data" / "raw_data" / "fails.fa"

    model.number(seqs)

    return model


@pytest.mark.parametrize(
    "scheme",
    ["imgt", "kabat"],
)
def test_files_are_identical(anarcii_model, tmp_path, scheme, pytestconfig):
    suffix = f"_{scheme}"

    expected_file = (
        pytestconfig.rootpath
        / "tests"
        / "data"
        / "expected_data"
        / f"fails{suffix}_expected.msgpack"
    )

    gen_object = from_msgpack_map(expected_file)
    expected = next(gen_object)

    test = anarcii_model.to_scheme(scheme)

    # Ensure both lists have the same length
    assert len(expected) == len(test), (
        f"Expected list length {len(expected)} but got {len(expected)}"
    )

    # Iterate over both lists concurrently
    for nm, expected_item, test_item in zip(
        list(expected.keys()), list(expected.values()), list(test.values())
    ):
        expected_number, expected_data = expected_item["numbering"], expected_item
        test_number, test_data = test_item["numbering"], test_item

        if expected_number and test_number:
            expected_number = [((x[0][0], x[0][1]), x[1]) for x in expected_number]

        else:
            # The number format is empty tuple for expected, and empty list for test
            expected_number = None
            test_number = None

        assert expected_number == test_number, (
            f"Numbering for {nm} is different! "
            f"Expected: {expected_number}, Got: {test_number}"
        )
        reference = pytest.approx(expected_data["score"], abs=0.5)
        assert test_data["score"] == reference, (
            f"Scores differ more than 0.5 for {nm}! "
            f"Expected: {expected_data['score']}, Got: {test_data['score']}"
        )
