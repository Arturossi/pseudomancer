import argparse
from .__init__ import __version__
from .pipeline import check_dependencies
from .pipeline import run_pseudomancer_pipeline
import os
import sys


def get_args():
    description = (
        "pseudomancer: a command line tool for reconstructing"
        + " pseudogenes in prokaryotic genomes using mmseqs2"
    )
    main_parser = argparse.ArgumentParser(
        description=description,
        prog="pseudomancer",
    )

    # i/o args
    io_opts = main_parser.add_argument_group("Input and output")
    io_opts.add_argument(
        "-g",
        "--genome",
        dest="genome_file",
        help=(
            "FASTA-format file containing genome of interest which will be "
            + " queried for potential pseudogenes"
        ),
        type=str,
        required=True,
        default=None,
    )
    io_opts.add_argument(
        "--genus",
        dest="genus",
        help=(
            "Genus name for downloading reference protein sequences from NCBI RefSeq "
            + "(e.g., 'Mycobacterium'). Proteins will be downloaded, clustered at 99%% identity, "
            + "and used as queries for homolog searches."
        ),
        type=str,
        required=True,
        default=None,
    )
    io_opts.add_argument(
        "-o",
        "--out_dir",
        dest="output_dir",
        help="directory for output files",
        type=str,
        required=True,
        default=None,
    )

    # search args
    search_opts = main_parser.add_argument_group("Search arguments")
    search_opts.add_argument(
        "-e",
        "--evalue",
        dest="e_value",
        help="e-value threshold for identifying homologs using mmseqs2 (default: 1e-5)",
        type=float,
        required=False,
        default=1e-5,
    )

    # main parser args
    main_parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    args = main_parser.parse_args()
    main_parser.set_defaults(func=run_pseudomancer_pipeline)

    return args


def main():
    args = get_args()

    check_dependencies()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)
    
    # Run the pipeline
    run_pseudomancer_pipeline(
        genus=args.genus,
        genome_file=args.genome_file,
        output_dir=args.output_dir,
        evalue=args.e_value
    )


if __name__ == "__main__":
    main()
