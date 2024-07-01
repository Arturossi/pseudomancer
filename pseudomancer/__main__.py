import argparse
from .__init__ import __version__
from .tblastn import check_tblastn_version
from .tblastn import run_tblastn
import os
import sys


class UltimateHelpFormatter(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    pass


def get_args():
    description = (
        "pseudomancer: a command line tool for reconstructing"
        + " pseudogenes in prokaryotic genomes"
    )
    main_parser = argparse.ArgumentParser(
        description=description,
        prog="pseudomancer",
        formatter_class=UltimateHelpFormatter,
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
        "-p",
        "--proteins",
        dest="proteins_file",
        help=(
            "FASTA-format file containing containing protein sequences to be"
            + " used as queries for homolog searches against the genome of interest"
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

    # tblastn args
    blast_opts = main_parser.add_argument_group("tblastn arguments")
    blast_opts.add_argument(
        "-e",
        "--evalue",
        dest="e_value",
        help="e-value threshold for identifying homologs using tblastn",
        type=float,
        required=True,
        default=1e-5,
    )

    # main parser args
    main_parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + __version__
    )

    args = main_parser.parse_args()
    main_parser.set_defaults(func=run_tblastn)

    return args


def main():
    args = get_args(sys.argv[1:])

    check_blast_version()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)


if __name__ == "__main__":
    main()
