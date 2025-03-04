from .utils import ATOM_RECORDS, THREE_TO_ONE, find_repeated_sequences


def renumber_pdb_with_anarcii(
    file_path,
    inner_seq_type: str = "antibody",
    inner_mode: str = "accuracy",
    inner_batch_size: int = 4,
    inner_cpu: bool = False,
):
    """
    Renumber PDB files

    Args:
    - file_path (str): Directory containing input PDB files.
    """
    from anarcii import Anarcii

    # The sequences will be passed as a list to the inner model.
    # This means that unknown mode can be run successfully and Classifii is bypassed
    # at the earlier stage in the pipeline __init__.py code.

    inner_model = Anarcii(
        seq_type=inner_seq_type,
        mode=inner_mode,
        batch_size=inner_batch_size,
        cpu=inner_cpu,
        verbose=False,
    )

    out_name = file_path.replace(".pdb", "_anarcii.pdb")
    try:
        sequences = {}
        with open(file_path) as pdb_file:
            for line in pdb_file:
                if line[:6].strip() in ATOM_RECORDS:
                    chain_id = line[21:22].strip()
                    residue_name = line[17:20].strip()
                    residue_number = int(line[22:26].strip())
                    insertion_code = line[26:27].strip()
                    residue_key = (chain_id, residue_number, insertion_code)

                    if chain_id not in sequences:
                        sequences[chain_id] = []

                    if (
                        not sequences[chain_id]
                        or sequences[chain_id][-1][0] != residue_key
                    ):
                        sequences[chain_id].append(
                            (residue_key, THREE_TO_ONE.get(residue_name, "X"))
                        )

        renumbering_scheme = {}  # Extract the new numbering scheme from ANARCII

        num = 1
        for chain_id, sequence_info in sequences.items():
            num += 1

            sequence = "".join([res[1] for res in sequence_info])

            # # check for repeats and take the first 190
            if find_repeated_sequences(sequence):
                print("Repeated seq of 190 residues found.")
                # TODO: Use actual repeated subseq and offset, not arbitrary 190.
                sequence = sequence[:190]

            seq = [(chain_id, sequence)]
            result = inner_model.number(seq)
            nums, alignment = result[0]
            chain_call = alignment["chain_type"]
            score = alignment["score"]

            seen_lines = set()
            if chain_call in [
                "H",
                "L",
                "K",  # AB
                "A",
                "B",
                "D",
                "G",
            ]:  # TCR
                if score < 19.0 and not (
                    (("23", " "), "C") in nums
                    and (("41", " "), "W") in nums
                    and (("104", " "), "C") in nums
                ):
                    print(
                        "PDB chain: ",
                        chain_id,
                        "Failed, low score and missing conserved residues. Score: ",
                        score,
                    )
                    print(len(sequence), sequence)
                    continue

                elif score < 19.0 and (
                    (("23", " "), "C") in nums
                    and (("41", " "), "W") in nums
                    and (("104", " "), "C") in nums
                ):
                    print(
                        "PDB chain: ",
                        chain_id,
                        "ANARCII Chain (Score): ",
                        chain_call,
                        f"({score})",
                        " Low score with conserved - check the sequence!",
                    )
                    print(len(sequence), sequence)

                else:
                    print(
                        "PDB chain: ",
                        chain_id,
                        "ANARCII Chain (Score): ",
                        chain_call,
                        f"({score})",
                    )

                nums = [x for x in nums if x[1] != "-"]

                anarcii_seq = "".join([x[1].strip() for x in nums])
                i2 = sequence.find(anarcii_seq)
                transform = int(nums[0][0][0]) - i2

                i3 = 1
                i4 = 1
                last = int(nums[-1][0][0])

                for i in range(len(sequence)):
                    current_line = sequence_info[i][0]
                    if current_line in seen_lines:
                        break
                    else:
                        seen_lines.add(current_line)
                        if i <= i2:
                            renumbering_scheme[sequence_info[i][0]] = (
                                i + transform,
                                " ",
                            )
                        elif i3 < len(nums):
                            renumbering_scheme[sequence_info[i][0]] = nums[i3][0]
                            i2 += 1
                            i3 += 1
                        else:
                            renumbering_scheme[sequence_info[i][0]] = (last + i4, " ")
                            i4 += 1

        with open(file_path) as infile, open(out_name, "w") as outfile:
            for line in infile:
                if line[:6].strip() in ATOM_RECORDS:
                    chain_id = line[21:22].strip()
                    original_residue_number = int(line[22:26].strip())
                    original_insertion_code = line[26:27].strip()
                    key = (chain_id, original_residue_number, original_insertion_code)

                    if key in renumbering_scheme:
                        new_residue_number, new_insertion_code = renumbering_scheme[key]
                        line = (
                            line[:22]
                            + f"{new_residue_number:>4}"
                            + f"{new_insertion_code:>1}"
                            + line[27:]
                        )

                outfile.write(line)

    except OSError as e:
        print(f"Failed to process {file_path}: {e}")
        print(f"Error occurred in processing chain: {chain_id}")
