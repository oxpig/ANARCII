# COPY ANARCII TOKENISER
from collections.abc import Iterable

import numpy as np
from numpy.typing import NDArray

non_standard_aa = set("BOJUZ")

mhc_nums = [str(x) for x in range(1, 251)] + [str(x) for x in range(1001, 1093)]


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

        # MHCs ======================================================
        elif self.vocab_type == "protein_mhc":
            self.vocab = [
                self.pad,
                self.start,
                self.end,
                self.skip,
                *([x.upper() for x in "acdefghiklmnpqrstvwXy"]),
            ]

        elif self.vocab_type == "number_mhc":
            self.vocab = [
                self.pad,
                self.start,
                self.end,
                self.skip,
                *mhc_nums,
                "X",
                *[
                    "B2M",
                    "CD1A",
                    "CD1B",
                    "CD1C",
                    "CD1D",
                    "CD1E",
                    "H2-AA",
                    "H2-AB",
                    "H2-D1",
                    "H2-D4",
                    "H2-DMA",
                    "H2-DOA",
                    "H2-DOB",
                    "H2-EA",
                    "H2-EB",
                    "H2-K",
                    "H2-L",
                    "H2-M",
                    "H2-Q",
                    "H2-T",
                    "HFE",
                    "HLA-A",
                    "HLA-B",
                    "HLA-C",
                    "HLA-DMA",
                    "HLA-DMB",
                    "HLA-DOA",
                    "HLA-DOB",
                    "HLA-DPA1",
                    "HLA-DPB1",
                    "HLA-DQA1",
                    "HLA-DQA2",
                    "HLA-DQB1",
                    "HLA-DQB2",
                    "HLA-DRA",
                    "HLA-DRB1",
                    "HLA-DRB3",
                    "HLA-DRB4",
                    "HLA-DRB5",
                    "HLA-E",
                    "HLA-F",
                    "HLA-G",
                    "MIC",
                    "MR1",
                    "TAP1",
                    "TAP2",
                ],  # Unpack the MHC labels
            ]

        else:
            raise ValueError(f"Vocab type {vocab_type} not supported")

        super().__init__()
