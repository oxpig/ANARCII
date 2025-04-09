from collections.abc import Iterable

import numpy as np
from numpy.typing import NDArray

non_standard_aa = set("BOJUZ")


class Tokeniser:
    def __init__(self):
        vocab = getattr(self, "vocab", [])
        self.tokens = np.array(vocab, dtype=object)
        self.char_to_int = {c: i for i, c in enumerate(vocab)}
        if "X" in vocab:
            for char in non_standard_aa:
                self.char_to_int[char] = self.char_to_int["X"]

    def encode(self, sequence: Iterable[str]) -> NDArray[np.int32]:
        # Replace non-standard amino acids with 'X'
        standardised_sequence: list[int] = [self.char_to_int[char] for char in sequence]
        return np.array(standardised_sequence, np.int32)


class NumberingTokeniser(Tokeniser):
    def __init__(self, vocab_type="protein"):
        self.vocab_type = vocab_type
        self.pad = "<PAD>"
        self.start = "<SOS>"
        self.end = "<EOS>"
        self.skip = "<SKIP>"

        # Antibodies ==================================================
        if self.vocab_type == "protein_antibody":
            self.vocab = [
                self.pad,
                self.start,
                self.end,
                self.skip,
                *([x.upper() for x in "acdefghiklmnpqrstvwXy"]),
            ]

        elif self.vocab_type == "number_antibody":
            self.vocab = [
                self.pad,
                self.start,
                self.end,
                self.skip,
                *list(range(1, 129)),
                "X",
                "H",
                "L",
                "K",
            ]

        # TCRs ======================================================
        elif self.vocab_type == "protein_tcr":
            self.vocab = [
                self.pad,
                self.start,
                self.end,
                self.skip,
                *([x.upper() for x in "acdefghiklmnpqrstvwXy"]),
            ]

        elif self.vocab_type == "number_tcr":
            self.vocab = [
                self.pad,
                self.start,
                self.end,
                self.skip,
                *list(range(1, 129)),
                "X",
                "A",
                "B",
                "G",
                "D",
            ]

        else:
            raise ValueError(f"Vocab type {vocab_type} not supported")

        super().__init__()
