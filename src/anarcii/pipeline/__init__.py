from __future__ import annotations

import sys
import time

from anarcii.classifii import Classifii
from anarcii.inference.model_runner import ModelRunner
from anarcii.inference.window_selector import WindowFinder
from anarcii.input_data_processing import Input, coerce_input, split_sequences
from anarcii.input_data_processing.sequences import SequenceProcessor
from anarcii.output_data_processing.convert_to_legacy_format import convert_output
from anarcii.output_data_processing.schemes import convert_number_scheme
from anarcii.pipeline.configuration import configure_cpus, configure_device

# Processing output
from anarcii.pipeline.methods import (
    print_initial_configuration,
    to_csv,
    to_dict,
    to_imgt_regions,
    to_json,
    to_text,
)

if sys.version_info >= (3, 12):
    from itertools import batched
else:
    from itertools import islice

    def batched(iterable, n, *, strict=False):
        # batched('ABCDEFG', 3) → ABC DEF G
        if n < 1:
            raise ValueError("n must be at least one")
        iterator = iter(iterable)
        while batch := tuple(islice(iterator, n)):
            if strict and len(batch) != n:
                raise ValueError("batched(): incomplete batch")
            yield batch


class Anarcii:
    """
    This class instantiates the models based on user input.

    Then it runs the number method, detecting input type.

    Number method does:
        * Checking of input sequence/file type.
        * Based on input it formats to a dict of {name:seq } - SequenceProcessor
        * Processed seqs are passed to model which uses ModelRunner class to perform
        autogressive inference steps.
        * Numbered seqs can be returned as a list, as well as be written to:
             csv,json, txt

    IF:
        * Very long list of seqs, or a long fasta file - the process is broken up
        into chunks and the outputs written to a text file in the working dir.

        * PDB file - detected and renumbered in-situ, returning file_anarcii.pdb

        * UNKNOWN model - a classifer model Classifii is called on partially processed
        input seqs. This detects if they are TCRs or Antibodies. Then runs the relevant
        model - returning the mixed list of both types.

    """

    def __init__(
        self,
        seq_type: str = "antibody",
        mode: str = "accuracy",
        batch_size: int = 32,
        cpu: bool = False,
        ncpu: int = -1,
        output_format: str = "simple",  # legacy for old ANARCI
        verbose: bool = False,
        max_seqs_len=1024 * 100,
    ):
        # need checks that all adhere before running code.
        self.seq_type = seq_type.lower()
        self.mode = mode.lower()
        self.batch_size = batch_size
        self.verbose = verbose
        self.cpu = cpu
        self.max_seqs_len = max_seqs_len

        self.output_format = output_format.lower()

        self._last_numbered_output = None
        # Has a conversion to a new number scheme occured?
        self._last_converted_output = None
        self._alt_scheme = None

        # Attach methods
        self.print_initial_configuration = print_initial_configuration.__get__(self)
        self.to_text = to_text.__get__(self)
        self.to_csv = to_csv.__get__(self)
        self.to_json = to_json.__get__(self)
        self.to_dict = to_dict.__get__(self)
        self.to_imgt_regions = to_imgt_regions.__get__(self)

        # Get device and ncpu config
        self.ncpu = configure_cpus(ncpu)
        self.device = configure_device(self.cpu, self.ncpu)
        self.print_initial_configuration()

    def number(self, seqs: Input):
        seqs: dict[str, str] = split_sequences(coerce_input(seqs), self.verbose)
        n_seqs = len(seqs)

        # If there is more than one chunk, we will need to serialise the output.
        serialise = n_seqs > self.max_seqs_len

        if self.verbose:
            print(f"Length of sequence list: {n_seqs}")
            n_chunks = n_seqs // self.max_seqs_len + 1
            print(
                f"Processing sequences in {n_chunks} chunks of {self.max_seqs_len} "
                "sequences."
            )
            begin = time.time()

        if self.seq_type == "unknown":
            classifii_seqs = Classifii(batch_size=self.batch_size, device=self.device)

        for i, chunk in enumerate(batched(seqs.items(), self.max_seqs_len), 1):
            chunk = dict(chunk)

            if self.verbose:
                print(f"Processing chunk {i} of {n_chunks}.")

            if self.seq_type == "unknown":
                # Classify the sequences as TCRs or antibodies.
                classified = classifii_seqs(chunk)

                if self.verbose:
                    n_antibodies = len(classified["antibody"])
                    n_tcrs = len(classified["tcr"])
                    print("### Ran antibody/TCR classifier. ###\n")
                    print(f"Found {n_antibodies} antibodies and {n_tcrs} TCRs.")

                # Combine the numbered sequences.
                numbered = {}
                for seq_type, sequences in classified.items():
                    numbered.update(self.number_with_type(sequences, seq_type))

            else:
                numbered = self.number_with_type(chunk, self.seq_type)

            # Restore the original input order to the numbered sequences.
            numbered = {key: numbered[key] for key in chunk}

            if serialise:
                # TODO: Finalise the serialisation implementation.
                ...
            else:
                self._last_numbered_output = numbered

        if self.verbose:
            end = time.time()
            runtime = round((end - begin) / 60, 2)
            print(f"Numbered {n_seqs} seqs in {runtime} mins. \n")

        return convert_output(
            ls=self._last_numbered_output,
            format=self.output_format,
            verbose=self.verbose,
        )

    def to_scheme(self, scheme="imgt"):
        # Check if there's output to save
        if self._last_numbered_output is None:
            raise ValueError("No output to convert. Run the model first.")

        else:
            converted_seqs = convert_number_scheme(self._last_numbered_output, scheme)
            print(f"Last output converted to {scheme}")

            # The problem is we cannot write over last numbered output
            # Instead, the converted scheme is written to a new object
            # This allows it to be written to json/text/csv
            self._last_converted_output = converted_seqs
            self._alt_scheme = scheme

            return converted_seqs

    def number_with_type(self, seqs: dict[str, str], seq_type):
        model = ModelRunner(
            seq_type, self.mode, self.batch_size, self.device, self.verbose
        )
        window_model = WindowFinder(seq_type, self.mode, self.batch_size, self.device)

        processor = SequenceProcessor(seqs, model, window_model, self.verbose)
        tokenised_seqs, offsets = processor.process_sequences()

        # Perform numbering.
        return model(tokenised_seqs, offsets)
