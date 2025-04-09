import re

import torch

from anarcii.inference.model_runner import ModelRunner
from anarcii.inference.window_selector import WindowFinder
from anarcii.input_data_processing import TokenisedSequence
from anarcii.input_data_processing.tokeniser import Tokeniser

from .utils import pick_windows, split_seq

# A regex pattern to match no more than 200 residues, containing a 'CWC' pattern
# (cysteine followed by 5–25 residues followed by a tryptophan followed by 50–80
# residues followed by another cysteine) starting no later than the 41st residue. The
# pattern greedily captures 0–40 residues (labelled 'start') preceding the CWC pattern,
# then gredily captures the CWC pattern (labelled 'cwc') in a lookahead.  The next
# string of up to 160 residues (labelled 'end') is also greedily captured in a
# lookahead.  The search poition is then advanced to just before the trailing C of the
# captured CWC pattern, effectively making the 'C...W...' search atomic.  This allows
# matches to overlap, except for the 'C...W...' sections of the CWC groups.  The desired
# string of up to 200 residues must be reconstructed by combining the 'start' and 'end'
# groups.
cwc_pattern = re.compile(
    r"""
        (?P<start>.{,40})                # Capture up to 40 residues.
        (?=(?P<cwc>                      # Zero-width search capturing a CWC pattern.
            (?P<atom>C.{5,25}W.{50,80})C # Prepare an atomic match to 'C...W...' of CWC.
        ))
        (?=(?P<end>.{,160}))             # Zero-width search capturing up to 160 chars.
        (?P=atom)                        # Move to the terminating C of the matched CWC.
    """,
    re.VERBOSE,
)


class SequenceProcessor:
    """
    This class takes a dict of sequences  {name: seq}. As well as pre-defined models
    that relate to the sequence type (antibody, TCR, shark).

    It has several steps it performs to pre-process the list of seqs so it can be
    consumed by the language model. These include:

    # 1
    * Checking for long seqs that exceed the context window (200 residues)
    * Working out what "window" within the long seq should be passed to the model.
    * holding the offsets to allow us to translate the indices back to the original
      long seq.

    # 2
    * Sorting the tuple by length of seqs to ensure we can pad batches of seqs that all
    share a similar length - to reduce unnecessary autoregressive infercence steps.

    # 3
    * Tokenising the sequences to numbers - then putting these into torch tensors.

    """

    def __init__(
        self,
        seqs: dict[str, str],
        model: ModelRunner,
        window_model: WindowFinder,
        verbose: bool,
    ):
        """
        Args:
            seqs (dict): A dictionary, keys are sequence IDs and values are sequences.
            model (torch.nn.Module): PyTorch model for processing full sequences.
            window_model (torch.nn.Module): modification of the above model that uses
            a one step decoder to get get a single logit value representing
            score for the input window (sequence fragment).
            verbose (bool): Whether to print detailed logs.
        """
        self.seqs: dict[str, str] = seqs
        self.model: ModelRunner = model
        self.window_model: WindowFinder = window_model
        self.verbose: bool = verbose
        self.offsets: dict[str, int] = {}

    def process_sequences(self):
        # Step 1: Handle long sequences
        self._handle_long_sequences()

        # Step 2: Sort sequences by length
        self._sort_sequences_by_length()

        # Step 3: Tokenize sequences
        return self._tokenize_sequences(), self.offsets

    def _handle_long_sequences(self):
        n_jump = 3
        long_seqs = {key: seq for key, seq in self.seqs.items() if len(seq) > 200}

        if long_seqs and self.verbose:
            print(
                f"\n {len(long_seqs)} Long sequences detected - running in sliding "
                "window. This is slow."
            )

        for key, sequence in long_seqs.items():
            # first try a simple regex to look for cwc
            cwc_matches = list(cwc_pattern.finditer(sequence))
            seq_strings = [m.group("start") + m.group("end") for m in cwc_matches]
            cwc_strings = [m.group("cwc") for m in cwc_matches]

            if cwc_matches:
                # Output the integer index of a high scoring window
                cwc_winner = pick_windows(cwc_strings, self.window_model)

                if cwc_winner is not None:
                    # Append the start offset
                    self.offsets[key] = cwc_matches[cwc_winner].start()
                    # Replace the input sequence
                    self.seqs[key] = seq_strings[cwc_winner]
                    # print(seq_strings[cwc_winner])
                    continue

            # No CWC match found proceed to window
            # If no cwc pattern is found, use the sliding window approach.
            # Split the sequence into 90-residue chunks and pick the best.
            windows = split_seq(sequence, n_jump=n_jump)

            best_window = pick_windows(windows, model=self.window_model, fallback=True)

            # Ensures start_index is at least 0.
            start_index = max((best_window * n_jump) - 40, 0)
            end_index = (best_window * n_jump) + 160

            # Append the start offset
            self.offsets[key] = start_index
            # Replace the input sequence
            self.seqs[key] = sequence[start_index:end_index]

        if long_seqs and self.verbose:
            print("Max probability windows selected.\n")

    def _sort_sequences_by_length(self):
        self.seqs = dict(sorted(self.seqs.items(), key=lambda x: len(x[1])))

    def _tokenize_sequences(self) -> dict[str, TokenisedSequence]:
        aa: Tokeniser = self.model.sequence_tokeniser
        tokenized_seqs = {}

        for name, seq in self.seqs.items():
            bookend_seq = [aa.start, *seq, aa.end]
            try:
                tokenized_seqs[name] = torch.from_numpy(aa.encode(bookend_seq))
            except KeyError as e:
                print(
                    f"Sequence could not be numbered. Contains an invalid residue: {e}"
                )
                tokenized_seqs[name] = torch.from_numpy(aa.encode(["F"]))

        return tokenized_seqs
