from __future__ import annotations

import sys
import time
import uuid
from itertools import chain, count, repeat
from pathlib import Path

import gemmi
import msgpack

from anarcii.classifii import Classifii
from anarcii.inference.model_runner import CUTOFF_SCORE, ModelRunner
from anarcii.inference.window_selector import WindowFinder
from anarcii.input_data_processing import Input, coerce_input, split_sequences
from anarcii.input_data_processing.sequences import SequenceProcessor
from anarcii.output_data_processing.convert_to_legacy_format import legacy_output
from anarcii.output_data_processing.schemes import convert_number_scheme
from anarcii.pipeline.configuration import configure_cpus, configure_device

# Processing output
from anarcii.pipeline.methods import (
    print_initial_configuration,
    to_csv,
    to_imgt_regions,
    to_json,
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


# Record all groups in mmCIF output, including _atom_site.auth_atom_id and .auth_comp_id
mmcif_output_groups = gemmi.MmcifOutputGroups(True, auth_all=True)
# We don't modify ATOM serial numbers and CONECT records.  Preserve them in PDB output.
pdb_write_options = gemmi.PdbWriteOptions(preserve_serial=True, conect_records=True)


def format_timediff(timediff: int | float) -> str:
    """
    Format a time difference in seconds as hours, minutes and seconds strings.

    Args:
        runtime:  The time difference in seconds.

    Returns:
        Time difference formatted as 'H hours, MM minutes, SS.SS seconds'.
    """
    hours, remainder = divmod(timediff, 3600)
    mins, secs = divmod(remainder, 60)

    hours = f"{hours:.0f} hr, " if hours else ""
    mins = f"{mins:{'02' if hours else ''}.0f} min, " if hours or mins else ""
    secs = f"{secs:{'02' if hours or mins else ''}.2f} sec"

    return f"{hours}{mins}{secs}"


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
        legacy_format: bool = False,  # legacy for old ANARCI
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

        self._last_numbered_output: dict | Path | None = None
        # Has a conversion to a new number scheme occured?
        self._last_converted_output = None
        self._alt_scheme = None

        # Attach methods
        self.print_initial_configuration = print_initial_configuration.__get__(self)
        self.to_csv = to_csv.__get__(self)
        self.to_json = to_json.__get__(self)
        self.to_imgt_regions = to_imgt_regions.__get__(self)

        # Get device and ncpu config
        self.ncpu = configure_cpus(ncpu)
        self.device = configure_device(self.cpu, self.ncpu)
        self.print_initial_configuration()

    def number(self, seqs: Input, legacy_format=False):
        seqs, structure = coerce_input(seqs)
        if not structure:
            # Do not split sequences on delimiter characters if the input was in PDBx or
            # PDB format.  We assume that PDBx/PDB files will have chains identified
            # individually.
            seqs: dict[str, str] = split_sequences(seqs, self.verbose)
        n_seqs = len(seqs)

        if self.verbose:
            print(f"Length of sequence list: {n_seqs}")
            n_chunks = -(n_seqs // -self.max_seqs_len)

            print(
                f"Processing sequences in {n_chunks} chunks of {self.max_seqs_len} "
                "sequences."
            )
            begin = time.time()

        if self.seq_type == "unknown":
            classifii_seqs = Classifii(batch_size=self.batch_size, device=self.device)

        # If there is more than one chunk, we will need to serialise the output.
        if serialise := n_seqs > self.max_seqs_len:
            id = uuid.uuid4()
            self._last_numbered_output = Path(f"anarcii-{id}-imgt.msgpack")

            # If we serialise we always need to tell the user.
            print(
                "\n",
                f"Serialising output to {self._last_numbered_output} as the number of "
                f"sequences exceeds the serialisation limit of {self.max_seqs_len}.\n",
            )

            packer = msgpack.Packer()
            # Initialise a MessagePack map with the expected number of sequences, so we
            # can later stream the key value pairs, rather than needing to create a
            # separate MessagePack map for each chunk.
            with self._last_numbered_output.open("wb") as f:
                f.write(packer.pack_map_header(n_seqs))

        for i, chunk in enumerate(batched(seqs.items(), self.max_seqs_len), 1):
            chunk = dict(chunk)

            if self.verbose:
                print(f"Processing chunk {i} of {n_chunks}.")

            if self.seq_type == "unknown":
                # Classify the sequences as TCRs or antibodies.
                classified = classifii_seqs(chunk)

                if self.verbose:
                    n_antibodies = len(classified.get("antibody", ()))
                    n_tcrs = len(classified.get("tcr", ()))
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

            # If the sequences came from a PDB(x) file, renumber them in the associated
            # data structure.
            if structure:
                for (model_index, chain_id), numbering in numbered.items():
                    if self.verbose:
                        print(f"PDBx model index, chain ID: {model_index}, {chain_id}")
                    if numbered_sequence_qa(numbering, self.verbose):
                        renumber_pdbx(structure, model_index, chain_id, numbering)

            if serialise:
                # Stream the key-value pairs of the results dict to the previously
                # initialised MessagePack map.
                with self._last_numbered_output.open("ab") as f:
                    for item in chain.from_iterable(numbered.items()):
                        f.write(packer.pack(item))
            else:
                self._last_numbered_output = numbered

        if self.verbose:
            end = time.time()
            print(f"Numbered {n_seqs} seqs in {format_timediff(end - begin)}.\n")

        # If our sequences came from a PDBx or PDB file, write a renumbered version.
        if structure:
            write_pdbx_file(structure)

        if legacy_format and not serialise:
            return legacy_output(self._last_numbered_output, verbose=self.verbose)
        else:
            return self._last_numbered_output

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


def numbered_sequence_qa(numbered: dict, verbose=False) -> bool:
    """
    Quality assurance check for an ANARCII-numbered sequence.

    Given an ANARCII-numbered sequence, check whether it is 'good' according to the
    following criteria:
      1. The assigned chain type is either 'HLK' (antibody) or 'ABDG' (TCR); and
      2. Either of the following criteria are met:
         a) The numbering model score is at least 19.
         b) The sequence contains conserved residues at the following IMGT positions:
            - 23: C
            - 41: W
            - 104: C

    Args:
        numbered:  An ANARCII-numbered seuence and associated metadata.

    Returns:
        bool:  True if the criteria are met.
    """
    if numbered["chain_type"] in (
        "HLK"  # Antibody
        "ABDG"  # TCR
    ):
        if verbose:
            print(
                f"ANARCII chain type (score): {numbered['chain_type']} "
                f"({numbered['score']})\n",
                f"Sequence length: {len(numbered['numbering'])}\n",
                f"Sequence: {numbered['numbering']}",
            )
        if numbered["score"] >= CUTOFF_SCORE:
            return True
        else:
            conserved_residues = {
                (("23", " "), "C"),
                (("41", " "), "W"),
                (("104", " "), "C"),
            }
            if conserved_residues.intersection(numbered["numbering"]):
                if verbose:
                    print("Low score with conserved residues — check the sequence!")
                return True
            else:
                return False
    else:
        return False


def renumber_pdbx(
    structure: gemmi.Structure, model_index: int, chain_id: str, numbered: dict
) -> None:
    """
    Write residue numbers from an ANARCII-numbered sequence to a Gemmi structure.

    Args:
        structure:    Representation of a PDBx or PDB file.
        model_index:  Index of the relevant model in the file.
        chain_id:     ID of the relevant chain in the model.
        numbering:    ANARCII model output for a given sequence.
    """
    # Get the sequence indicated by the model index and chain ID.
    polymer: gemmi.ResidueSpan = structure[model_index][chain_id].get_polymer()
    # Drop gap marks ('-') from the numbered sequence.  They do not exist in the file.
    no_gaps = ((num, res) for num, res in numbered["numbering"] if res != "-")
    # Get the residue numbering and one-letter peptide sequence as separate tuples.
    numbers, sequence = zip(*no_gaps)
    # Find the number of the first numbered residue.
    (first_number, _), *_, (last_number, _) = numbers

    try:
        # Get the numbering offset, by matching the numbered sequence to the original...
        offset: int = polymer.make_one_letter_sequence().index("".join(sequence))
    except ValueError:
        # ... or by falling back on the model's reported start index.
        offset: int = numbered["query_start"]

    # Generate numbers for the residues in the file that precede the numbered sequence.
    backward_fill = zip(range(first_number - offset, first_number), repeat(" "))
    forward_fill = zip(count(last_number + 1), repeat(" "))
    numbers = chain(backward_fill, numbers, forward_fill)

    # Residue by residue, write the new numbering.
    for residue, number in zip(structure[model_index][chain_id], numbers):
        residue.seqid = gemmi.SeqId(*number)


def write_pdbx_file(structure: gemmi.Structure, scheme="imgt") -> None:
    """
    Write a Gemmi PDBx structure to file.

    Use the same format as the source file, as determined by `structure.input_format`.
    Label the file with `structure.name` and the name of the numbering scheme used.

    Args:
        structure:  Representation of a PDBx or PDB file.
        scheme:     Numbering scheme used to generate the structure.
    """
    stem = f"{structure.name.lower()}-anarcii-{scheme}"

    if structure.input_format is gemmi.CoorFormat.Pdb:
        structure.write_pdb(f"{stem}.pdb", pdb_write_options)

    else:
        document = structure.make_mmcif_document(mmcif_output_groups)

        if structure.input_format is gemmi.CoorFormat.Mmcif:
            document.write_file(f"{stem}.cif")

        elif structure.input_format is gemmi.CoorFormat.Mmjson:
            with open(f"{stem}.json", "w") as f:
                f.write(document.as_json(mmjson=True))
