import sys
from pathlib import Path

import gemmi
import pytest

from anarcii import Anarcii

if sys.version_info < (3, 11):
    import os
    from contextlib import contextmanager

    @contextmanager
    def chdir(path):
        """
        Change directory for the duration of the context.

        A non-reentrant poor man's `contextlib.chdir`, for Python versions < 3.11.
        """
        _old_cwd = os.getcwd()
        try:
            yield os.chdir(path)
        finally:
            os.chdir(_old_cwd)
else:
    from contextlib import chdir


raw_filenames = "1kb5.pdb", "8kdm.pdb"
raw_paths = [Path(f) for f in raw_filenames]
reference_paths = [p.with_stem(f"{p.stem}_anarcii") for p in raw_paths]
test_paths = [p.with_stem(f"{p.stem}-anarcii-imgt") for p in raw_paths]

model = Anarcii(
    seq_type="unknown",
    batch_size=1,
    cpu=True,
    ncpu=16,
    mode="accuracy",
    verbose=True,
)


@pytest.fixture(scope="session")
def anarcii_model(tmp_path_factory, pytestconfig) -> Path:
    """
    Renumber some source PDB riles and return the path to their temporary directory.
    """
    tmp_path = tmp_path_factory.mktemp("renumbered-pdbs-")
    raw_data = pytestconfig.rootpath / "tests" / "data" / "raw_data"
    with chdir(tmp_path):
        for filename in raw_paths:
            # At present, PDB renumbering only writes output PDBx and PDB files to the
            # working directory.
            model.number(raw_data / filename)

    return tmp_path


@pytest.mark.parametrize(
    "reference,test", zip(reference_paths, test_paths), ids=raw_filenames
)
def test_files_are_identical(pytestconfig, anarcii_model, reference, test):
    """Generate and check renumbered PDB files."""
    expected_data = pytestconfig.rootpath / "tests" / "data" / "expected_data"

    expected_structure = gemmi.read_structure(str(expected_data / reference))
    expected_structure.setup_entities()
    test_structure = gemmi.read_structure(str(anarcii_model / test))
    test_structure.setup_entities()

    expected_numbering = [
        [[residue.seqid for residue in chain] for chain in model]
        for model in expected_structure
    ]
    test_numbering = [
        [[residue.seqid for residue in chain] for chain in model]
        for model in test_structure
    ]
    assert test_numbering == expected_numbering
