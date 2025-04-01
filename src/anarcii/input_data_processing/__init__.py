from __future__ import annotations

import gzip
import re
from collections.abc import Iterator
from functools import partial
from itertools import chain
from pathlib import Path
from typing import TypeAlias

import gemmi
import torch

# Valid user input types.
Input: TypeAlias = (
    Path | str | tuple[str, str] | list[str | tuple[str, str]] | dict[str, str]
)
# A TokenisedSequence is a torch.Tensor of dtype np.int32.
TokenisedSequence: TypeAlias = torch.Tensor


gz_suffixes = {".gz", ".z"}
# Supported FASTA file suffixes.  Peptide sequences only, no nucleotides.
fasta_suffixes = {".fasta", ".fas", ".fsa", ".fa", ".faa", ".mpfa"}
# Supported PIR file suffixes.
pir_suffixes = {".pir", ".nbrf", ".ali"}

supported_extensions = fasta_suffixes | pir_suffixes

paired_sequence_delimiters = r"-\/"
split_pattern = paired_sequence_delimiters.replace("\\", r"\\")
split_pattern = re.compile(rf"[{split_pattern}]")


def file_input(path: Path) -> dict[str, str]:
    """
    Extract peptide sequence strings from a file.

    Supported file formats are:
    * FASTA (.fasta, .fas, .fa, .faa, .mpfa and their gzipped equivalents).
    * NBRF/PIR (.pir, .nbrf, .ali and their gzipped equivalents).

    Args:
        path (pathlib.Path): Path to the input file.

    Returns:
        dict[str, str]: The values are sequence strings and the keys are
                        names/descriptions, which are assumed to be unique in the input
                        file.
    """
    if (fasta_suffixes | pir_suffixes).intersection(path.suffixes):
        with gzip.open(path, "rt") if path.suffix in gz_suffixes else open(path) as f:
            entries: list[gemmi.FastaSeq] = gemmi.read_pir_or_fasta(f.read())

        return {e.header: e.seq for e in entries if e.header and e.seq}

    else:
        raise ValueError(
            f"{path.name} has an unsupported file extension.  These are supported:\n"
            f"{', '.join(sorted(supported_extensions))}\n"
            "and gzipped equivalents (*.gz, *.z)."
        )


def coerce_input(input_data: Input) -> dict[str, str]:
    """
    Coerce varied input sequence data formats into a dictionary.

    Accepts one or more peptide sequence strings, packaged up in a variety of ways,
    producing a dictionary with the sequences as values and the accompanying labels as
    keys:
    * A file path, which will be read to extract the sequence strings and their labels;
    * A single sequence string, which will be labelled with the key `sequence`;
    * A list of sequence strings which will be labelled sequentially with keys
      `sequence-1`, `sequence-2`, etc.;
    * A dictionary of name-sequence pairs, which will be returned unmodified;
    * A tuple of name-sequence pairs, or a list thereof, will be converted into a
      dictionary of the same.

    Args:
        input_data (Input): Peptide sequence strings, optionally labelled with names.

    Raises:
        TypeError: An unrecognised type of input data was provided.

    Returns:
        dict[str, str]: The values are sequence strings and the keys are names, which
                        are generated if not provided.
    """
    try:
        # Capture the cases list[tuple[str, str]] | dict[str, str],
        # containing name-sequence pairs.
        return dict(input_data)

    # The only non-iterable sub-type of Input is pathlib.Path.
    except TypeError:
        # Capture the case of file input (pathlib.Path).
        if isinstance(input_data, Path):
            return file_input(input_data)

    # The remaining sub-types of Input are str | tuple[str, str] | list[str].
    except ValueError:
        if isinstance(input_data, str):
            # Capture the case of file input (str).
            path = Path(input_data)
            if path.suffix:
                return file_input(path)

            # Capture the case of a single peptide sequence (str).
            return {"Sequence": input_data}

        if isinstance(input_data, tuple):
            # Capture the case of a single name-sequence pair (tuple[str, str]).
            name, sequence = input_data
            return {name: sequence}

        if isinstance(input_data, list):
            # Capture the case of a list of peptide sequences (list[str]), labelling
            # sequentially with 'Sequence 1', 'Sequence 2', etc..
            width = len(str(len(input_data)))
            return {
                f"Sequence {i:0{width}d}": seq for i, seq in enumerate(input_data, 1)
            }

    raise TypeError("Invalid input type.")


def split_sequence(
    name: str, sequence: str, verbose: bool = False
) -> Iterator[tuple[str, str]]:
    """
    Split a sequence string on any of several standard delimiter characters.

    Sequence delimiter characters are stripped from the start and end of a sequence
    string.

    Args:
        name:      Sequence name.
        sequence:  Peptide sequence, possibly containing delimiters.
        verbose:   If true, print a warning message when the sequence is split.

    Returns:
        For split sequences, an iterator of 2-tuple name-sequence pairs.  For an unsplit
        sequence, a 1-tuple containing a 2-tuple name-sequence pair.
    """
    # Strip leading and trailing delimiters.
    sequence = sequence.strip(paired_sequence_delimiters)
    # Check for remaining delimiters.
    if "-" in sequence or "\\" in sequence or "/" in sequence:
        if verbose:
            delimiters_str = "' or '".join(paired_sequence_delimiters)
            print(
                f"'{delimiters_str}' found in sequence {name}.",
                "Assuming this is a paired sequence and splitting into parts.",
            )

        # Split the sequence on these delimiters.
        split_parts = re.split(split_pattern, sequence)
        width = len(str(len(split_parts)))
        # Create named parts
        for i, part in enumerate(split_parts, 1):
            yield f"{name}-{i:0{width}d}", part
    else:
        # If no delimiters, yield the sequence as is.
        yield name, sequence


def split_sequences(seqs: dict[str, str], verbose: bool = False) -> dict[str, str]:
    """
    Split sequence strings on any of several standard delimiter characters.

    Sequence delimiter characters are stripped from the start and end of each sequence
    string.

    Args:
        seqs:     Dictionary of name-sequence pairs.
        verbose:  If true, print a warning message for each split sequence.

    Returns:
        A dictionary in which sequence strings that have been split on delimiters have
        been replaced with the split parts, each labelled with a unique name.  The input
        key order is retained in the output, aside from the split sequences, which are
        inserted in order at the position of their originating input sequence.
    """
    splitter = partial(split_sequence, verbose=verbose)
    return dict(chain.from_iterable(map(splitter, seqs.keys(), seqs.values())))
