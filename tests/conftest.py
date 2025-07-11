import os
import pytest
import tempfile
import shutil
from pathlib import Path


@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_genome_file(temp_dir):
    """Create a sample genome FASTA file for testing."""
    genome_file = os.path.join(temp_dir, "test_genome.fasta")
    with open(genome_file, 'w') as f:
        f.write(">test_contig_1\n")
        f.write("ATGAAACGTACGCTGAAATTTGCGCAGTAG\n")
        f.write(">test_contig_2\n") 
        f.write("ATGGCAACCGTAGGCTGATCGATCGACGACTAG\n")
    return genome_file


@pytest.fixture
def sample_protein_file(temp_dir):
    """Create a sample protein FASTA file for testing."""
    protein_file = os.path.join(temp_dir, "test_proteins.faa")
    with open(protein_file, 'w') as f:
        f.write(">test_protein_1\n")
        f.write("MKRLLAISLLLAVVTSLLAAPYVKA\n")
        f.write(">test_protein_2\n")
        f.write("MATAIGDRSTLTA\n")
    return protein_file