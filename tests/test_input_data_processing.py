import gzip

import pytest

from anarcii.input_data_processing import (
    coerce_input,
    fasta_suffixes,
    file_input,
    gz_suffixes,
    pir_suffixes,
)


@pytest.fixture(scope="class")
def empty_input_files(tmp_path_factory):
    """Create empty input files, of all valid types."""
    tmp_path = tmp_path_factory.mktemp("file_types")
    fasta_paths = [
        tmp_path / f"foo{suffix}" for suffix in fasta_suffixes | pir_suffixes
    ]
    for path in fasta_paths:
        path.write_text(">")

    fasta_gz_paths = (
        path.with_suffix(f"{path.suffix}{gz_suffix or ''}")
        for path in fasta_paths
        for gz_suffix in gz_suffixes
    )
    for path in fasta_gz_paths:
        with gzip.open(path, "w") as f:
            f.write(b">")

    return tmp_path


cases = {
    "single-sequence-string": (
        "test",
        {"Sequence": "test"},
    ),
    "single-name-sequence-tuple": (
        ("name", "test"),
        {"name": "test"},
    ),
    "list-single-sequence-string": (
        ["test"],
        {"Sequence 1": "test"},
    ),
    "list-multiple-sequence-strings": (
        ["test-1", "test-2"],
        {"Sequence 1": "test-1", "Sequence 2": "test-2"},
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


@pytest.fixture
def sample_fasta(request):
    return request.config.rootpath / "tests" / "data" / "raw_data" / "100_seqs.fa"


class TestFiles:
    @pytest.mark.parametrize("compression", {None} | gz_suffixes)
    @pytest.mark.parametrize("suffix", fasta_suffixes | pir_suffixes)
    def test_valid_empties(self, empty_input_files, suffix, compression):
        """Check that correctly formatted empty files are read as empty dictionaries."""
        filename = f"foo{suffix}{compression or ''}"
        assert file_input(empty_input_files / filename) == ({}, None)

    def test_invalid_extension(self, empty_input_files):
        """Check that an invalid file extension raises a ValueError."""
        with pytest.raises(ValueError):
            file_input(empty_input_files / "foo.bar")


@pytest.mark.parametrize("input_data, expected", cases.values(), ids=cases)
def test_coerce_input(input_data, expected):
    """Check that non-file input data is coerced into the expected format."""
    assert coerce_input(input_data) == (expected, None)


def test_invalid_input_type():
    """Check that an invalid input type raises a TypeError."""
    with pytest.raises(TypeError) as exc_info:
        coerce_input(0)

    assert str(exc_info.value) == "Invalid input type."


@pytest.mark.parametrize("input_function", [file_input, coerce_input])
def test_file_input(sample_fasta, input_function):
    """Check that the first sequecence in a valid FASTA file is read correctly."""
    seqs, structure = input_function(sample_fasta)
    assert seqs["sequence0"] == (
        "SETLSLTCSVYGASISNSNSYWGWIRQPPGKRLEWLGSIYDSGSTSYNPSLSS"
        "RVTISVDTSKNQVSLRLNSVTAADTGVYYCARHRDPPGSRWIFYYYYMDLWG"
    )
    assert structure is None


def test_file_input_str_or_path(sample_fasta):
    """Check that the file is read equally well when a file name string is passed."""
    assert coerce_input(sample_fasta) == coerce_input(str(sample_fasta))
