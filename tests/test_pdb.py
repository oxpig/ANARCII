import filecmp
import pathlib
import shutil
import sys

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
raw_file_paths = map(pathlib.Path, raw_filenames)
output_filenames = [p.with_stem(f"{p.stem}_anarcii") for p in raw_file_paths]

model = Anarcii(
    seq_type="unknown",
    batch_size=1,
    cpu=False,
    ncpu=16,
    mode="accuracy",
    verbose=True,
)


@pytest.fixture(scope="session")
def anarcii_model(tmp_path_factory, pytestconfig) -> pathlib.Path:
    """
    Renumber some source PDB riles and return the path to their temporary directory.
    """
    tmp_path = tmp_path_factory.mktemp("renumbered-pdbs-")
    raw_data = pytestconfig.rootpath / "tests" / "data" / "raw_data"
    with chdir(tmp_path):
        for filename in raw_filenames:
            # At present, PDB renumbering only works in place, acting on a file in the
            # working directory.  Accordingly, copy the source data into our temporary
            # directory.
            shutil.copy2(raw_data / filename, tmp_path)
            model.number(filename)

    return tmp_path


@pytest.mark.parametrize("filename", output_filenames)
def test_files_are_identical(pytestconfig, anarcii_model, filename):
    """Generate and check both text and json files."""
    expected_data = pytestconfig.rootpath / "tests" / "data" / "expected_data"
    expected_file = expected_data / filename
    test_file = anarcii_model / filename

    assert filecmp.cmp(expected_file, test_file, shallow=False), (
        f"Files {expected_file} and {test_file} are different!"
    )
