__all__ = ["ATOM_RECORDS", "THREE_TO_ONE", "find_repeated_sequences", "repeat_pattern"]

# Constants
import re

ATOM_RECORDS = "ATOM"
THREE_TO_ONE = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}

# Minimum sequence length for detecting a repeat.
MIN_SEQUENCE_LENGTH = 50

# Pattern for detecting repeated sequences of 50 or more residues.
repeat_pattern = re.compile(rf"(.{{{MIN_SEQUENCE_LENGTH},}})(?=.*?\1)")


def find_repeated_sequences(sequence):
    """
    Counts the number of repeated subsequences of at least 50 residues in a sequence.

    Parameters:
    - sequence: The sequence string to search in.

    Returns:
    - An integer representing the number of repeated subsequences.
    """

    return len(repeat_pattern.findall(sequence))
