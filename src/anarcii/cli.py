import argparse
import sys

from anarcii.pipeline import Anarcii

parser = argparse.ArgumentParser(
    description="Run the Anarcii model on sequences or a fasta file."
)
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


def main():
    # Define the argument parser
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

    try:
        out = model.number(args.input)
    except TypeError as e:
        sys.exit(str(e))

    if not args.output:
        for name, query in out.items():
            # Print to screen
            print(
                f" ID: {name}\n",
                f"Chain: {query['chain_type']}\n",
                f"Score: {query['score']}\n",
                f"Error: {query['error']}",
            )
            print({"".join(map(str, n)).strip(): res for n, res in query["numbering"]})

    elif args.output.endswith(".csv"):
        model.to_csv(args.output)
    elif args.output.endswith(".json"):
        model.to_json(args.output)
    # TODO:  Add msgpack support.
    # elif args.output.endswith(".msgpack"):
    #     model.to_msgpack(args.output)
    else:
        raise ValueError("Output file must end in .csv, or .json.")


if __name__ == "__main__":
    main()
