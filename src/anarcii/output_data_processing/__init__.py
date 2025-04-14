from __future__ import annotations

import csv
from collections.abc import Iterable, Iterator
from itertools import chain, pairwise
from pathlib import Path
from typing import BinaryIO, TextIO, TypeAlias

from sortedcontainers import SortedSet

from anarcii.utils import _open_msgpack_map_file

NumberedResidue: TypeAlias = tuple[tuple[int, str], str]
NumberedResidues: TypeAlias = list[NumberedResidue] | tuple[NumberedResidue, ...]

# For IMGT, insertions are numbered in reverse lexicographic order at these positions.
imgt_reversed = 33, 61, 112


# Minimal CSV columns.
metadata_columns = "Name", "Chain", "Score", "Query start", "Query end"
required_residue_numbers = {(n, " ") for n in range(1, 129)}


def numbered_sequence_dict(numbering: NumberedResidues) -> dict[str, str]:
    """
    Convert a list or tuple of numbered residues to a dictionary.

    Each numbering `tuple[int, str]` of residue number and insertion character is
    coerced into a string key by concatenating the integer and string parts, stripping
    blank insertion characters.  The residue letter is taken as the corresponding value.

    Args:
        numbering: A list or tuple of tuples of the form
                   ((residue number, insertion character), residue letter)

    Returns:
        A dictionary with the concatenated numbering strings as keys and the residue
        letters as values
    """
    return {str(num) + ins.strip(): res for (num, ins), res in numbering}


def _imgt_order_segments(
    numbers: Iterable[tuple[int, str]],
) -> Iterator[Iterator[tuple[int, str]]]:
    """
    Sort IMGT residue numbers, taking into account reversed numbering for insertions.

    Args:
        numbers: A SortedSet of residue number strings.

    Yields:
        An iterable of ordered residue number strings.
    """
    numbers = SortedSet(numbers)
    half_open = True, False
    for low, high in pairwise((None, *imgt_reversed)):
        # The first range is open-ended, so we use None as the lower bound.
        yield numbers.irange((low + 1,) if low else None, (high,), inclusive=half_open)
        # Reverse the insertions in the latter half of each CDR.
        yield numbers.irange((high,), (high + 1,), inclusive=half_open, reverse=True)

    # Finally, yield all the remaining numbers after the CDR3 insertion region.
    yield numbers.irange(minimum=(high + 1,))


def imgt_order(numbers: Iterable[tuple[int, str]]) -> Iterator[tuple[int, str]]:
    """
    Sort IMGT residue numbers, taking into account reversed numbering for insertions.

    Args:
        numbers: A SortedSet of residue number strings.

    Returns:
        An iterable of ordered residue number strings.
    """
    # Sort the segments and concatenate them
    return chain.from_iterable(_imgt_order_segments(numbers))


def write_csv(numbered: dict, path: Path | str) -> None:
    """
    Write an ANARCII model results dictionary to a CSV file.

    The results dictionary may contain multiple numbered sequences, which will be
    aligned when written to the CSV file.  The file will contain the following columns:
    - Name: The name of the sequence.
    - Chain: The sequence's chain type ('F' in the case of a failure).
    - Score: The model's score for its numbering of the sequence.
    - Query start: The position of the first residue numbered by the model (this is left
      blank if the model failed to number the sequence).
    - Query end: The position of the last residue numbered by the model (this is left
      blank if the model failed to number the sequence).
    - One column for each residue number present in any of the sequences.  Residue
      numbers 1–128 are always included.  Residue columns are sorted in ascending number
      order, except in the case of IMGT numbering, where the system of inward numbering
      of CDR insertions is respected.

    In the table of sequences, residues are represented by their one-letter codes, or by
    '-' for absences.

    Args:
        numbered:  An ANARCII model results dictionary.
        path:      The path at which to write the CSV file.
    """
    residue_numbers = required_residue_numbers

    rows = []
    for name, result in numbered.items():
        numbering = result["numbering"] or []
        residue_numbers.update(number for number, _ in numbering)

        rows.append(
            {
                "Name": name,
                "Chain": result["chain_type"],
                "Score": result["score"],
                "Query start": result["query_start"],
                "Query end": result["query_end"],
                **numbered_sequence_dict(numbering),
            }
        )

    # Assume all sequences use the same scheme.  In any case, there's no point aligning
    # multiple sequences if they have have been numbered using different schemes.
    if result["scheme"] == "imgt":
        # Reverse certain insertions as necessary for IMGT numbering.
        residue_numbers = imgt_order(residue_numbers)
    else:
        residue_numbers = sorted(residue_numbers)

    residue_columns = (str(num) + ins.strip() for num, ins in residue_numbers)
    columns = [*metadata_columns, *residue_columns]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, restval="-")
        writer.writeheader()
        writer.writerows(rows)


def _stream_msgpack_file_to_csv_file(
    f: BinaryIO, g: TextIO, chunk_size: int = 100 * 1024
) -> None:
    """
    Stream serialised ANARCII model results from MessagePack to CSV.

    Standard ANARCII model results are read from a MessagePack file containing a single
    map.  The results map may contain multiple numbered sequences, which will be aligned
    when written to the CSV file.  Sequences will be streamed from the MessagePack map
    to the CSV file in batches of `chunk_size` sequences.

    The CSV file will contain the following columns:
    - Name: The name of the sequence.
    - Chain: The sequence's chain type ('F' in the case of a failure).
    - Score: The model's score for its numbering of the sequence.
    - Query start: The position of the first residue numbered by the model (this is left
      blank if the model failed to number the sequence).
    - Query end: The position of the last residue numbered by the model (this is left
      blank if the model failed to number the sequence).
    - One column for each residue number present in any of the sequences.  Residue
      numbers 1–128 are always included.  Residue columns are sorted in ascending number
      order, except in the case of IMGT numbering, where the system of inward numbering
      of CDR insertions is respected.

    In the table of sequences, residues are represented by their one-letter codes, or by
    '-' for absences.

    Args:
        f:           A file object for the MessagePack input.  Must be opened in binary
                     mode and contain a MessagePack map as the first entry.
        g:           A file object for the output.  Must be opened in text mode with
                     `newline=''`.
        chunk_size:  Streaming chunk size.  Number of sequences to read, convert and
                     write at a time.
    """
    residue_numbers = required_residue_numbers

    # A first pass over the MessagePack map to collect all residue numbers.
    for results in _open_msgpack_map_file(f, chunk_size):
        for result in results.values():
            residue_numbers.update(number for number, _ in (result["numbering"] or []))

    # Assume all sequences use the same scheme.  In any case, there's no point aligning
    # multiple sequences if they have have been numbered using different schemes.
    if result["scheme"] == "imgt":
        # Reverse certain insertions as necessary for IMGT numbering.
        residue_numbers = imgt_order(residue_numbers)
    else:
        residue_numbers = sorted(residue_numbers)

    residue_columns = (str(num) + ins.strip() for num, ins in residue_numbers)
    columns = [*metadata_columns, *residue_columns]

    writer = csv.DictWriter(g, fieldnames=columns, restval="-")
    writer.writeheader()

    # A second pass over the input iterable to write the numbered sequences to the file.
    f.seek(0)
    for results in _open_msgpack_map_file(f, chunk_size):
        rows = [
            {
                "Name": name,
                "Chain": result["chain_type"],
                "Score": result["score"],
                "Query start": result["query_start"],
                "Query end": result["query_end"],
                **numbered_sequence_dict(result["numbering"] or []),
            }
            for name, result in results.items()
        ]

        writer.writerows(rows)


def stream_msgpack_to_csv(
    msgpack_map_path: Path | str, csv_path: Path | str, chunk_size: int = 100 * 1024
) -> None:
    """
    Stream serialised ANARCII model results from MessagePack to CSV.

    Standard ANARCII model results are read from a MessagePack file containing a single
    map.  The results map may contain multiple numbered sequences, which will be aligned
    when written to the CSV file.  Sequences will be streamed from the MessagePack map
    to the CSV file in batches of `chunk_size` sequences.

    The CSV file will contain the following columns:
    - Name: The name of the sequence.
    - Chain: The sequence's chain type ('F' in the case of a failure).
    - Score: The model's score for its numbering of the sequence.
    - Query start: The position of the first residue numbered by the model (this is left
      blank if the model failed to number the sequence).
    - Query end: The position of the last residue numbered by the model (this is left
      blank if the model failed to number the sequence).
    - One column for each residue number present in any of the sequences.  Residue
      numbers 1–128 are always included.  Residue columns are sorted in ascending number
      order, except in the case of IMGT numbering, where the system of inward numbering
      of CDR insertions is respected.

    In the table of sequences, residues are represented by their one-letter codes, or by
    '-' for absences.

    Args:
        msgpack_map_path:  Path to a MessagePack file.  It must contain a map as the
                           first entry, representing an ANARCII model result dictionary.
        path:              Path for the output CSV file.
        chunk_size:        Streaming chunk size.  Number of sequences to read, convert
                           and write at a time.
    """
    with open(msgpack_map_path, "rb") as f, open(csv_path, "w", newline="") as g:
        _stream_msgpack_file_to_csv_file(f, g, chunk_size)
