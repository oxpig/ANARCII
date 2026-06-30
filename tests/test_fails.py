import pytest

from anarcii.output_data_processing.schemes import convert_number_scheme


@pytest.mark.parametrize(
    "scheme",
    ["imgt", "kabat"],
)
def test_failed_sequence_is_preserved_without_mutating_source(scheme):
    failure = {
        "numbering": None,
        "chain_type": "F",
        "score": 0.0,
        "scheme": "imgt",
        "query_start": None,
        "query_end": None,
        "error": "Could not apply numbering",
    }
    numbered_sequences = {"failed_sequence": failure.copy()}

    converted = convert_number_scheme(numbered_sequences, scheme)

    assert converted == {
        "failed_sequence": {
            **failure,
            "scheme": scheme,
        }
    }
    assert numbered_sequences == {"failed_sequence": failure}
    assert converted["failed_sequence"] is not numbered_sequences["failed_sequence"]
