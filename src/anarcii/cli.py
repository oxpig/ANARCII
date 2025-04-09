import argparse

from anarcii.pipeline import Anarcii


def main():
    # Define the argument parser
    parser = argparse.ArgumentParser(
        description="Run the Anarcii model on sequences or a fasta file."
    )

    # Add command-line flags and options
    parser.add_argument(
        "input", type=str, help="Input sequence as a string or path to a fasta file."
    )
    parser.add_argument(
        "-t",
        "--seq_type",
        type=str,
        default="antibody",
        help="Sequence type (default: antibody).",
    )
    parser.add_argument(
        "-b",
        "--batch_size",
        type=int,
        default=512,
        help="Batch size for processing (default: 512).",
    )
    parser.add_argument(
        "-c", "--cpu", action="store_true", help="Run on CPU (default: False)."
    )
    parser.add_argument(
        "-n",
        "--ncpu",
        type=int,
        default=-1,
        help="Number of CPU threads to use (default: 1).",
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        default="accuracy",
        choices=["accuracy", "speed"],
        help="Mode for running the model (default: accuracy).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Specify the output file (must end in .txt, .csv or .json).",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose output."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Initialize the model
    model = Anarcii(
        seq_type=args.seq_type,
        batch_size=args.batch_size,
        cpu=args.cpu,
        ncpu=args.ncpu,
        mode=args.mode,
        verbose=args.verbose,
    )

    # Check if the input is a file or a single sequence
    if args.input.endswith((".fasta", ".fa", ".fa.gz", ".fasta.gz")):
        print(f"Processing fasta file: {args.input}")
        out = model.number(args.input)
    elif args.input.endswith(".pdb"):
        print(f"Processing PDB file: {args.input}")
        out = model.number(args.input)
    else:
        print(f"Processing sequence: {args.input}")
        out = model.number([args.input])

    if not args.output:
        for i in range(len(out)):
            # Print to screen
            print(
                " ID: ",
                out[i][1]["query_name"],
                "\n",
                "Chain: ",
                out[i][1]["chain_type"],
                "\n",
                "Score: ",
                out[i][1]["score"],
                "\n",
                "Error: ",
                out[i][1]["error"],
            )
            [print(x) for x in out[i][0]]
    elif args.output.endswith(".csv"):
        model.to_csv(args.output)
    elif args.output.endswith(".txt"):
        model.to_txt(args.output)
    elif args.output.endswith(".json"):
        model.to_json(args.output)
    else:
        raise ValueError("Output file must end in .txt, .csv, or .json.")


if __name__ == "__main__":
    main()
