from __future__ import annotations

import csv
from collections.abc import Iterable, Iterator
from itertools import chain, pairwise, repeat
from pathlib import Path
from typing import TextIO, TypeAlias

from sortedcontainers import SortedSet

NumberedResidue: TypeAlias = tuple[tuple[int, str], str]
NumberedResidues: TypeAlias = list[NumberedResidue] | tuple[NumberedResidue, ...]

# For IMGT, insertions are numbered in reverse lexicographic order at these positions.
imgt_reversed = 33, 61, 112


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


def _imgt_order_segments(numbers: SortedSet) -> Iterator[Iterable[str]]:
    """
    Sort IMGT residue numbers, taking into account reversed numbering for insertions.

    Args:
        numbers: A SortedSet of residue number strings.

    Yields:
        An iterable of ordered residue number strings.
    """
    half_open = True, False
    for low, high in pairwise((None, *imgt_reversed)):
        # The first range is open-ended, so we use None as the lower bound.
        yield numbers.irange((low + 1,) if low else None, (high,), inclusive=half_open)
        # Reverse the insertions in the latter half of each CDR.
        yield numbers.irange((high,), (high + 1,), inclusive=half_open, reverse=True)

    # Finally, yield all the remaining numbers after the CDR3 insertion region.
    yield numbers.irange(minimum=(high + 1,))


def imgt_order(numbers: SortedSet) -> Iterable[str]:
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
    metadata_columns = "Name", "Chain", "Score", "Query start", "Query end"
    required_residue_numbers = zip(range(1, 129), repeat(" "))
    residue_numbers = SortedSet(required_residue_numbers)

    rows = []
    for name, result in numbered.items():
        numbering = result.get("numbering", [])
        residue_numbers.update(number for number, _ in numbering)

        rows.append(
            {
                "Name": name,
                "Chain": result["chain_type"],
                "Score": result["score"],
                "Query start": result.get("query_start"),
                "Query end": result.get("query_end"),
                **numbered_sequence_dict(numbering),
            }
        )

    # Assume all sequences use the same scheme.  In any case, there's no point aligning
    # multiple sequences if they have have been numbered using different schemes.
    if result["scheme"] == "imgt":
        # Reverse certain insertions as necessary for IMGT numbering.
        residue_numbers = imgt_order(residue_numbers)

    residue_columns = (str(num) + ins.strip() for num, ins in residue_numbers)
    columns = chain(metadata_columns, residue_columns)

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, restval="-")
        writer.writeheader()
        writer.writerows(rows)


def _stream_csv_to_file(numbered: Iterable[dict], f: TextIO) -> None:
    """
    Stream an iterable of ANARCII model results dictionary to a single CSV file.

    The results dictionaries may each contain multiple numbered sequences, which will be
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
        numbered:  An iterable of ANARCII model results dictionaries.
        f:         A file object for the output.  Must be opened in text mode with
                   `newline=''`.
    """
    metadata_columns = "Name", "Chain", "Score", "Query start", "Query end"
    required_residue_numbers = zip(range(1, 129), repeat(" "))
    residue_numbers = SortedSet(required_residue_numbers)

    # A first pass over the input iterable to collect all residue numbers.
    for results in numbered:
        for result in results.values():
            residue_numbers.update(number for number, _ in result.get("numbering", []))

    # Assume all sequences use the same scheme.  In any case, there's no point aligning
    # multiple sequences if they have have been numbered using different schemes.
    if result["scheme"] == "imgt":
        # Reverse certain insertions as necessary for IMGT numbering.
        residue_numbers = imgt_order(residue_numbers)

    residue_columns = (str(num) + ins.strip() for num, ins in residue_numbers)
    columns = chain(metadata_columns, residue_columns)

    writer = csv.DictWriter(f, fieldnames=columns, restval="-")
    writer.writeheader()

    # A second pass over the input iterable to write the numbered sequences to the file.
    for results in numbered:
        rows = []
        for name, result in results.items():
            numbering = result.get("numbering", [])
            residue_numbers.update(number for number, _ in numbering)

            rows.append(
                {
                    "Name": name,
                    "Chain": result["chain_type"],
                    "Score": result["score"],
                    "Query start": result.get("query_start"),
                    "Query end": result.get("query_end"),
                    **numbered_sequence_dict(numbering),
                }
            )

        writer.writerows(rows)


def stream_csv(numbered: Iterable[dict], path: Path | str) -> None:
    """
    Stream an iterable of ANARCII model results dictionary to a single CSV file.

    The results dictionaries may each contain multiple numbered sequences, which will be
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
        numbered:  An iterable of ANARCII model results dictionaries.
        path:      The path at which to write the CSV file.

    """
    with open(path, "w", newline="") as f:
        _stream_csv_to_file(numbered, f)
