import subprocess
import sys


def check_blast_version(blast_path: str = "tblastn"):
    """
    Check the version of BLAST+ installed on the system.

    Parameters
    ----------
    blast_path : str
        Path to the BLAST+ executable. If not provided, the function will
        search the system PATH for the executable.
    """
    try:
        blast_version = subprocess.check_output(
            [blast_path, "-version"],
            stderr=subprocess.STDOUT,
            shell=True,
        )
    except subprocess.CalledProcessError:
        sys.stderr.write(
            "Could not find BLAST+ executable. Please ensure that BLAST+ is"
            + " installed and in your PATH."
        )
        sys.exit(1)

    return blast_version
