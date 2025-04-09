import ast
import json

from anarcii.output_data_processing.list_to import (
    return_imgt_regions,
    write_csv,
    write_json,
)


def print_initial_configuration(self):
    """Print initial configuration details if verbose mode is enabled."""
    if self.verbose:
        print(f"Batch size: {self.batch_size}")
        print(
            "\tSpeed is a balance of batch size and length diversity. "
            "Adjust accordingly.\n",
            "\tSeqs all similar length (+/-5), increase batch size. "
            "Mixed lengths (+/-30), reduce.\n",
        )
        if not self.cpu:
            if self.batch_size < 512:
                print(
                    "Consider a batch size of at least 512 for optimal GPU performance."
                )
            elif self.batch_size > 512:
                print("For A100 GPUs, a batch size of 1024 is recommended.")
        else:
            print("Recommended batch size for CPU: 8.")


def to_csv(self, file_path):
    # Check if there's output to save
    if self._last_numbered_output is None:
        raise ValueError("No output to save. Run the model first.")

    elif self._last_converted_output and not self.max_len_exceed:
        write_csv(self._last_converted_output, file_path)
        print(
            f"Last output saved to {file_path} in alternate scheme: {self._alt_scheme}."
        )

    elif self._last_converted_output and self.max_len_exceed:
        raise ValueError(
            f"Cannot renumber more than {1024 * 100} sequences and convert"
            " to alternate scheme. Feature update coming soon!"
        )

    elif self.max_len_exceed:
        print(
            "Writing to a aligned CSV file may use a lot of RAM for millions of "
            "sequences, consider to_text(filepath) or to_json(filepath) for memory-"
            "efficient solutions."
        )
        with open(self.text_) as file:
            loaded_data = [ast.literal_eval(line.strip()) for line in file]

        write_csv(loaded_data, file_path)
        print(f"Last output saved to {file_path}")

    else:
        write_csv(self._last_numbered_output, file_path)
        print(f"Last output saved to {file_path}")


def to_json(self, file_path):
    # Check if there's output to save
    if self._last_numbered_output is None:
        raise ValueError("No output to save. Run the model first.")

    elif self._last_converted_output and not self.max_len_exceed:
        write_json(self._last_converted_output, file_path)
        print(
            f"Last output saved to {file_path} in alternate scheme: {self._alt_scheme}."
        )

    elif self._last_converted_output and self.max_len_exceed:
        raise ValueError(
            f"Cannot renumber more than {1024 * 100} sequences and convert"
            " to alternate scheme. Feature update coming soon!"
        )

    elif self.max_len_exceed:
        with open(self.text_) as infile, open(file_path, "w") as outfile:
            # Start the JSON array
            outfile.write("[\n")

            first = True  # To manage commas between JSON objects
            for line in infile:
                try:
                    # Parse each line safely
                    data = ast.literal_eval(line.strip())
                except (ValueError, SyntaxError):
                    print(f"Skipping invalid line: {line.strip()}")
                    continue

                # Write the parsed line as JSON, adding commas where needed
                if not first:
                    outfile.write(",\n")
                first = False
                json.dump(data, outfile)

            # End the JSON array
            outfile.write("\n]\n")

        print(f"Last output saved to {file_path}")

    else:
        write_json(self._last_numbered_output, file_path)
        print(f"Last output saved to {file_path}")


def to_imgt_regions(self):
    # Check if there's output to save
    if self._last_numbered_output is None:
        raise ValueError("No output. Run the model first.")

    elif self.max_len_exceed:
        with open(self.text_) as file:
            loaded_data = [ast.literal_eval(line.strip()) for line in file]
        ls = return_imgt_regions(loaded_data)
        return ls

    else:
        ls = return_imgt_regions(self._last_numbered_output)
        return ls
