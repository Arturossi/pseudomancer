import argparse
from .__init__ import __version__
from .tblastn import check_tblastn_version
from .tblastn import run_tblastn


class UltimateHelpFormatter(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    pass


def get_options():
    description = "pseudomancer: s command line tool for reconstructing "+
    " pseudogenes in prokaryotic genomes
    main_parser = argparse.ArgumentParser(
        description = description,
        prog = "pseudomancer",
        formatter_class=UltimateHelpFormatter
    )

    # i/o args
    io_opts = parser.add_argument_group('Input and output')
    io_opts.add_argument(
        "-g",
        "--genome",
        dest="genome_file",
        required=True,
        help=("input genome file containing putative pseudogenes"),
        type=str)
    io_opts.add_argument(
        "-p",
        "--proteins",
        dest="proteins_file",
        required=True,
        help=("fasta file containing containing protein sequences"),
        type=str)
    io_opts.add_argument(
        "-o",
        "--out_dir",
        dest="output_dir",
        required=True,
        help="directory for output files",
        type=str
    )

    # tblastn args
    tb_opts = parser.add_argument_group('tblastn arguments')
    tb_opts.add_argument(
        "-e,
        "--evalue",
        dest="e_value",
        required=True,
        help=("e-value threshold for identifying homologs using tblastn"),
        type=str)

    # main parser args
    main_parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + __version__
    )
    # parse arguments and run function
    args = main_parser.parse_args()
    args.func(args)

    return

def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo('MGEfinder version: ' + __version__)
    ctx.exit()

def main():
    args = get_options(sys.argv[1:])


if __name__ == "__main__":
    cli()
