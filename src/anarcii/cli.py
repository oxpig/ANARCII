import argparse
import sys

from anarcii import __version__
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
    choices=["antibody", "tcr", "vnar", "vhh", "shark", "unknown"],
    help=(
        "Sequence type to process: antibody, tcr, vnar/vhh/shark or unknown"
        "(default: antibody)."
    ),
)
parser.add_argument(
    "-b",
    "--batch_size",
    type=int,
    default=512,
    metavar="N",
    help="Batch size for processing (default: 512).",
)
parser.add_argument(
    "-c",
    "--cpu",
    action="store_true",
    help="Run on CPU only, even if a GPU is available.",
)
parser.add_argument(
    "-n",
    "--ncpu",
    type=int,
    default=-1,
    metavar="N",
    help="Number of CPU threads to use.  If -1 (the default), ANARCII will use one "
    "thread per available CPU core.",
)
parser.add_argument(
    "--max_seqs_len",
    type=int,
    default=102400,
    metavar="N",
    help=(
        "Maximum number of sequences to process before moving to batch mode and "
        "saving the numbered sequences in MessagePack file (default: 102400)."
    ),
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
    "--scheme",
    type=str,
    choices=["martin", "kabat", "chothia", "imgt", "aho"],
    default="imgt",
    help=(
        "Numbering scheme to use: martin, kabat, chothia, imgt, or aho (default: imgt)."
    ),
)
parser.add_argument(
    "-o",
    "--output",
    type=str,
    default=None,
    metavar="FILE",
    help="Specify the output file (must end in .csv or .msgpack).",
)
parser.add_argument(
    "-v", "--verbose", action="store_true", help="Enable verbose output."
)
parser.add_argument(
    "-V", "--version", action="version", version=f"%(prog)s {__version__}"
)


def main(args=None):
    args = parser.parse_args(args)

    # Initialize the model
    model = Anarcii(
        seq_type=args.seq_type,
        batch_size=args.batch_size,
        cpu=args.cpu,
        ncpu=args.ncpu,
        mode=args.mode,
        verbose=args.verbose,
        max_seqs_len=args.max_seqs_len,
    )

    try:
        model.number(args.input)
        out = model.to_scheme(args.scheme)
    except TypeError as e:
        sys.exit(str(e))

    if not args.output:
        for name, query in out.items():
            # Print to screen
            print(
                f" ID: {name}\n",
                f"Chain: {query['chain_type']}\n",
                f"Score: {query['score']}\n",
                f"Query start: {query['query_start']}\n",
                f"Query end: {query['query_end']}\n",
                f"Scheme: {query['scheme']}\n",
                f"Error: {query['error']}",
            )
            print({"".join(map(str, n)).strip(): res for n, res in query["numbering"]})

    elif args.output.endswith(".csv"):
        model.to_csv(args.output)
    elif args.output.endswith(".msgpack"):
        model.to_msgpack(args.output)
    else:
        raise ValueError("Output file must end in .csv, or .json.")


if __name__ == "__main__":
    main()
