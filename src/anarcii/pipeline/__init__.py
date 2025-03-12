import os
import time

from anarcii.classifii import Classifii
from anarcii.classifii.utils import join_mixed_types
from anarcii.inference.model_runner import ModelRunner
from anarcii.inference.window_selector import WindowFinder

# Classes
from anarcii.input_data_processing.sequences import SequenceProcessor
from anarcii.output_data_processing.convert_to_legacy_format import convert_output
from anarcii.output_data_processing.schemes import convert_number_scheme
from anarcii.pdb_process import renumber_pdb_with_anarcii
from anarcii.pipeline.batch_process import batch_process
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

# Functions for processing input
from anarcii.pipeline.utils import (
    count_lines_with_greater_than,
    is_tuple_list,
    read_fasta,
    split_sequence,
)


# This is the orchestrator of the whole pipeline.
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
        # For use with SCFVs or artificial constructs:
        scfv_or_concatenated_chains: bool = False,
        batch_size: int = 8,
        cpu: bool = False,
        ncpu: int = -1,
        output_format: str = "simple",  # legacy for old ANARCI
        verbose: bool = False,
        debug_input=False,
        max_seqs_len=1024 * 100,
    ):
        # need checks that all adhere before running code.
        self.seq_type = seq_type.lower()
        self.mode = mode.lower()
        self.batch_size = batch_size
        self.verbose = verbose
        self.cpu = cpu
        self.text_ = "output.txt"
        self.max_seqs_len = max_seqs_len
        self.max_len_exceed = False

        self.output_format = output_format.lower()
        self.scfv = scfv_or_concatenated_chains
        self.debug_input = debug_input
        self.unknown = False

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

        # shark model
        self.shark_model = ModelRunner(
            "shark", self.mode, self.batch_size, self.device, self.verbose
        )
        self.shark_window = WindowFinder(
            "shark", self.mode, self.batch_size, self.device, self.scfv
        )

        # Antibody model
        self.ig_model = ModelRunner(
            "antibody", self.mode, self.batch_size, self.device, self.verbose
        )
        self.ig_window = WindowFinder(
            "antibody", self.mode, self.batch_size, self.device, self.scfv
        )

        # TCR model
        self.tcr_model = ModelRunner(
            "tcr", self.mode, self.batch_size, self.device, self.verbose
        )
        self.tcr_window = WindowFinder(
            "tcr", self.mode, self.batch_size, self.device, self.scfv
        )

    def number(self, seqs, warnings=False):
        if self.seq_type.lower() == "unknown" and not (
            ".pdb" in seqs or ".mmcif" in seqs
        ):
            # classify the  sequences into tcrs or antibodies
            classifii_seqs = Classifii(batch_size=self.batch_size, device=self.device)

            if isinstance(seqs, list):
                if is_tuple_list(seqs):
                    dict_of_seqs = {}
                    for t in seqs:
                        split_seqs = split_sequence(t[0], t[1], self.verbose)
                        dict_of_seqs.update(split_seqs)

                elif not is_tuple_list(seqs):
                    dict_of_seqs = {}
                    for n, t in enumerate(seqs):
                        name = f"seq{n}"
                        # Split each sequence as needed and update the dictionary
                        split_seqs = split_sequence(name, t, self.verbose)
                        dict_of_seqs.update(split_seqs)

            elif isinstance(seqs, str) and os.path.exists(seqs) and ".fa" in seqs:
                fastas = read_fasta(seqs, self.verbose)
                dict_of_seqs = {t[0]: t[1] for t in fastas}

            else:
                print(
                    "Input is not a list of sequences, nor a valid path to a fasta file"
                    "(must end in .fa or .fasta), nor a pdb file (.pdb)."
                )
                return []

            list_of_tuples_pre_classifii = list(dict_of_seqs.items())
            list_of_names = list(dict_of_seqs.keys())
            antibodies, tcrs = classifii_seqs(list_of_tuples_pre_classifii)

            if self.verbose:
                print("### Ran antibody/TCR classifier. ###\n")
                print(f"Found {len(antibodies)} antibodies and {len(tcrs)} TCRs.")

            antis_out = self.number_with_type(
                antibodies, "antibody", chunk=False, warnings=warnings
            )
            # If max length has been exceeded here (but not in the next chunk).
            # You need to ensure that the TCR numberings are written to the same file.
            # check status of self.max_len_exceed.
            if self.max_len_exceed:
                chunk_subsequent = True
            else:
                chunk_subsequent = False

            # We need to stay in unknown mode and append to an output file.
            self.unknown = True

            tcrs_out = self.number_with_type(
                tcrs, "tcr", chunk=chunk_subsequent, warnings=warnings
            )
            self.unknown = False  # Reset to false.

            self._last_numbered_output = join_mixed_types(
                antis_out, tcrs_out, list_of_names
            )

            return convert_output(
                ls=self._last_numbered_output,
                format=self.output_format,
                verbose=self.verbose,
            )

        # PDB files can be run directly, bypassing classifii at this stage.
        # They will be classified into ab/tcrs by an inner model which takes a list of
        #  seqs read from a PDB file (these can now pass through the Classifii code).
        else:
            self._last_numbered_output = self.number_with_type(
                seqs, self.seq_type, chunk=False, warnings=warnings
            )
            return convert_output(
                ls=self._last_numbered_output,
                format=self.output_format,
                verbose=self.verbose,
            )

    def to_scheme(self, scheme="imgt"):
        # Check if there's output to save
        if self._last_numbered_output is None:
            raise ValueError("No output to convert. Run the model first.")

        elif self.max_len_exceed:
            raise ValueError(
                f"Cannot renumber more than {1024 * 100} sequences and convert"
                " to alternate scheme. Feature update coming soon!"
            )
        else:
            converted_seqs = convert_number_scheme(self._last_numbered_output, scheme)
            print(f"Last output converted to {scheme}")

            # The problem is we cannot write over last numbered output
            # Instead, the converted scheme is written to a new object
            # This allows it to be written to json/text/csv
            self._last_converted_output = converted_seqs
            self._alt_scheme = scheme

            return converted_seqs

    def number_with_type(self, seqs, inner_type, chunk=False, warnings=False):
        if inner_type == "shark":
            model = self.shark_model
            window_model = self.shark_window
        elif inner_type == "antibody":
            model = self.ig_model
            window_model = self.ig_window
        elif inner_type == "tcr":
            model = self.tcr_model
            window_model = self.ig_window
        else:
            print("Error in defining sequence type to number.")

        # Reset this
        self.max_len_exceed = False

        # clear the output file
        if hasattr(self, "text_") and os.path.exists(self.text_) and not self.unknown:
            os.remove(self.text_)

        chunk_list = []
        begin = time.time()

        if isinstance(seqs, list):
            if is_tuple_list(seqs):
                if self.verbose:
                    print(
                        "Running on a list of tuples of format: [(name,sequence), ...]."
                    )

                nms = [x[0] for x in seqs]
                if len(nms) != len(set(nms)):  # Detect duplicates
                    raise SystemExit(
                        "Error: Duplicate names found."
                        + "Please ensure all names are unique."
                    )

                dict_of_seqs = {}
                for t in seqs:
                    # Split each sequence as needed and update the dictionary
                    split_seqs = split_sequence(t[0], t[1], self.verbose)
                    dict_of_seqs.update(split_seqs)

            elif not is_tuple_list(seqs):
                if self.verbose:
                    print("Running on a list of strings: [sequence, ...].")

                dict_of_seqs = {}
                for n, t in enumerate(seqs):
                    name = f"seq{n}"
                    # Split each sequence as needed and update the dictionary
                    split_seqs = split_sequence(name, t, self.verbose)
                    dict_of_seqs.update(split_seqs)

            if self.verbose:
                print("Length of sequence list: ", len(dict_of_seqs))

            # If the list is huge - breakup into chunks of 1M.
            if len(dict_of_seqs) > self.max_seqs_len or chunk:
                print(
                    "\nMax # of seqs exceeded.",
                    f"Running chunks of {self.max_seqs_len}.\n",
                )

                keys = list(dict_of_seqs.keys())  # Convert dictionary keys to a list
                num_seqs = len(keys)

                num_chunks = (len(dict_of_seqs) // self.max_seqs_len) + 1
                for i in range(num_chunks):
                    chunk_keys = keys[
                        i * self.max_seqs_len : (i + 1) * self.max_seqs_len
                    ]
                    chunk_list.append({k: dict_of_seqs[k] for k in chunk_keys})

                self.max_len_exceed = True

        # Fasta file
        elif isinstance(seqs, str) and os.path.exists(seqs) and ".fa" in seqs:
            if self.verbose:
                print("Running on fasta file.")

            num_seqs = count_lines_with_greater_than(seqs)
            if self.verbose:
                print("Length of sequence list: ", num_seqs)

            if num_seqs > self.max_seqs_len:
                print(
                    f"Max # of seqs exceeded. Running chunks of {self.max_seqs_len}.\n"
                )
                num_chunks = (num_seqs // self.max_seqs_len) + 1
                for i in range(num_chunks):
                    fastas = read_fasta(seqs, self.verbose)[
                        i * self.max_seqs_len : (i + 1) * self.max_seqs_len
                    ]
                    chunk_list.append({t[0]: t[1] for t in fastas})
                self.max_len_exceed = True

            else:
                fastas = read_fasta(seqs, self.verbose)
                dict_of_seqs = {t[0]: t[1] for t in fastas}

        # PDB files
        elif (
            isinstance(seqs, str)
            and os.path.exists(seqs)
            and (".pdb" in seqs or ".mmcif" in seqs)
        ):
            # Unknown mode is taken care of here. Do not worry about passing classifii.
            print(f"Renumbering a PDB/mmCIF file in {self.seq_type} mode")
            numbered_chains = renumber_pdb_with_anarcii(
                seqs,
                inner_seq_type=self.seq_type,
                inner_mode=self.mode,
                inner_batch_size=self.batch_size,
                inner_cpu=self.cpu,
            )
            return numbered_chains

        else:
            print(
                "Input is not a list of sequences, nor a valid path to a fasta file "
                "(must end in .fa or .fasta), nor a pdb file (.pdb)."
            )
            return []

        if len(chunk_list) > 1:
            numbered_seqs = batch_process(
                chunk_list, model, window_model, self.verbose, self.text_, warnings
            )

            end = time.time()
            runtime = round((end - begin) / 60, 2)

            if self.verbose:
                print(f"Numbered {num_seqs} seqs in {runtime} mins. \n")

            print(
                f"\nOutput written to {self.text_}. Convert to csv or text with: "
                "model.to_csv(filepath) or model.to_text(filepath)"
            )

            return numbered_seqs

        else:
            # instantiate the Sequences class and process
            sequences = SequenceProcessor(
                dict_of_seqs, model, window_model, self.verbose, self.scfv
            )
            processed_seqs, offsets = sequences.process_sequences()

            # Offset for longseqs only - replace the indices...
            # ==========================================================================
            numbered_seqs = model(processed_seqs, offsets, warnings=warnings)
            # ==========================================================================

            end = time.time()
            runtime = round((end - begin) / 60, 2)

            if self.verbose:
                print(f"Numbered {len(numbered_seqs)} seqs in {runtime} mins. \n")

            if self.debug_input:
                print("Printing list of input seqs.")
                print(dict_of_seqs)

            return numbered_seqs
