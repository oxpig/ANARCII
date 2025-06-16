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


def padded_indices(count):
    width = len(str(count))
    return (f"{i + 1:0{width}d}" for i in range(count))


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
        scfv: bool,
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
            scfv (bool): whether to run in SCFV mode which looks for multiple IG/TCR
            regions in one sequence.
        """
        self.seqs: dict[str, str] = seqs
        self.model: ModelRunner = model
        self.window_model: WindowFinder = window_model
        self.verbose: bool = verbose
        self.offsets: dict[str, int] = {}
        self.scfv: bool = scfv

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

            if self.scfv:
                # Set up SCFV specific variables here. These can be played with.
                # The parameters below work best for most SCFV seqs tested.

                SCFV_JUMP = 1  # How many residues we increment along a sequence.
                SCFV_WINDOW_SIZE = 125  # Number of residues being scored.
                SCFV_WINDOW_NUM = int(SCFV_WINDOW_SIZE / SCFV_JUMP)

                SHIFT = int(50 / SCFV_JUMP)  # no of windows to move along: 50 residues
                SCFV_THRESHOLD = 20  # Score cut off for a given window

                windows = split_seq(
                    sequence, n_jump=SCFV_JUMP, window_size=SCFV_WINDOW_SIZE
                )

                data = pick_windows(
                    windows, model=self.window_model, scfv=True, fallback=True
                )

                ### Start by indentifying the minima - the sequence positions between
                # two regions which the model suggests contains IG/TCR content.
                minima = []
                start_idx = 0
                last_start = 0

                # Create a copy of data for later - we will reduce the size of data as
                # we iteratively search for the minima in the next 125 residues.
                probs = data

                # iterate through data and find minima that adhear to our conditions.
                while len(data) > 1:
                    min_value = min(data[:SCFV_WINDOW_NUM])
                    # The minima must be global...
                    # And not at the end of the sequence..
                    # Or the start...
                    if (
                        (min_value < SCFV_THRESHOLD)
                        and (
                            (len(probs) - (data.index(min_value) + last_start))
                            > 10 / SCFV_JUMP
                        )
                        and (data.index(min_value) + last_start) > 10 / SCFV_JUMP
                    ):
                        start_idx = data.index(min_value)
                        minima.append(start_idx + last_start)

                    # We need to ensure similar minima are not too close together
                    # move on 5 windows (SHIFT)
                    data = data[start_idx + SHIFT :]
                    last_start += start_idx + SHIFT

                    continue

                # Found windows >>> modify the seqs dict (still works for non SCFVs)
                # Remove the original key if it already exists
                seq = self.seqs[key]

                self.offsets.pop(key, None)
                self.seqs.pop(key, None)

                ### NOW LOOK FOR PEAKS (> threshold & within 50 residues of minima).
                offset = 0
                minima = minima + [len(probs)]

                idx = 1
                found = 0
                for i in range(len(minima)):
                    # For we want maxima in first 50 residues.
                    window = probs[offset : (offset + int(50 / SCFV_JUMP))]

                    # print("Offset:", offset,
                    #       "MAX:", max(window),
                    #       "Max IDX", probs.index(max(window)),
                    #       "Len:", len(window))

                    # Shift the offset to the new minima.
                    offset = minima[i]
                    if window and max(window) > SCFV_THRESHOLD:
                        # Add 2 to the peak index to ensure we contain the Ig domain
                        # This was simply found by trial and error (SORRY, no magic).
                        peak_idx_plus2 = probs.index(max(window)) + 2
                        new_key = f"{key}-{idx}"

                        # We will cap all sequences at 180.
                        if i == 0:
                            # First: must include sequence from 0 index.
                            window = seq[0 : (peak_idx_plus2 * SCFV_JUMP + 180)][:180]
                        else:
                            window = seq[
                                (peak_idx_plus2 * SCFV_JUMP) : (
                                    peak_idx_plus2 * SCFV_JUMP + 180
                                )  # increment by 180 aa
                            ]

                        self.offsets[new_key] = peak_idx_plus2 * SCFV_JUMP
                        self.seqs[new_key] = window
                        if self.verbose:
                            print(
                                f"Identified potential domain. Renamed to: {new_key}\n",
                                f"{window}",
                            )
                        idx += 1
                        found += 1
                if self.verbose:
                    print("")

                if found == 1:
                    print(f"Only found 1 domain for {key}.\n")
                elif found == 0:
                    print(f"Failed to find any domain for {key}.\n")
                elif found > 2:
                    print(f"Found more than 2 domains for {key}.\n")

                # Return the domains found to the user and do not enter the CWC loop.
                # CWC method was not robust to SCFVs during testing.
                continue

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
