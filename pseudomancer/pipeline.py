import subprocess
import sys
import os
import shutil
import tempfile
import requests
import pandas as pd
from pathlib import Path
from typing import List, Optional
import zipfile


def check_dependencies():
    """Check if required dependencies (mmseqs2, datasets) are available."""
    dependencies = ["mmseqs", "datasets"]
    missing = []
    
    for dep in dependencies:
        if not shutil.which(dep):
            missing.append(dep)
    
    if missing:
        sys.stderr.write(
            f"Missing required dependencies: {', '.join(missing)}\n"
            "Please install the following tools:\n"
            "- mmseqs2: https://github.com/soedinglab/MMseqs2\n"
            "- datasets: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/\n"
        )
        sys.exit(1)


def download_ncbi_assemblies(genus: str, output_dir: str) -> str:
    """
    Download NCBI RefSeq assemblies for a given genus.
    
    Parameters
    ----------
    genus : str
        Genus name (e.g., 'Mycobacterium')
    output_dir : str
        Directory to store downloaded assemblies
        
    Returns
    -------
    str
        Path to the merged protein FASTA file
    """
    assemblies_dir = os.path.join(output_dir, f"{genus.lower()}_assemblies")
    os.makedirs(assemblies_dir, exist_ok=True)
    
    # Use datasets to download assemblies
    print(f"Downloading {genus} assemblies from NCBI RefSeq...")
    cmd = [
        "datasets", "download", "genome", "taxon", genus,
        "--annotated", "--reference", "--assembly-level", "complete",
        "--include", "genome,protein,gff3,gtf,seq-report",
        "--filename", os.path.join(assemblies_dir, f"{genus.lower()}_dataset.zip")
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Failed to download assemblies: {e.stderr}\n")
        sys.exit(1)
    
    # Extract and merge protein files
    protein_files = []
    zip_path = os.path.join(assemblies_dir, f"{genus.lower()}_dataset.zip")
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(assemblies_dir)
    
    # Find all protein.faa files
    for root, dirs, files in os.walk(assemblies_dir):
        for file in files:
            if file == "protein.faa":
                protein_files.append(os.path.join(root, file))
    
    # Merge protein files
    merged_proteins = os.path.join(output_dir, f"{genus.lower()}_proteins.faa")
    with open(merged_proteins, 'w') as outfile:
        for protein_file in protein_files:
            with open(protein_file, 'r') as infile:
                outfile.write(infile.read())
    
    print(f"Merged {len(protein_files)} protein files into {merged_proteins}")
    return merged_proteins


def cluster_proteins_mmseqs2(protein_file: str, output_dir: str, identity: float = 0.99) -> str:
    """
    Cluster proteins using mmseqs2 at specified identity threshold.
    
    Parameters
    ----------
    protein_file : str
        Path to input protein FASTA file
    output_dir : str
        Directory for output files
    identity : float
        Sequence identity threshold (default: 0.99)
        
    Returns
    -------
    str
        Path to representative sequences FASTA file
    """
    base_name = os.path.splitext(os.path.basename(protein_file))[0]
    output_prefix = os.path.join(output_dir, f"{base_name}_clustered")
    
    print(f"Clustering proteins at {identity*100}% identity...")
    
    cmd = [
        "mmseqs", "easy-cluster",
        protein_file,
        output_prefix,
        output_dir,
        "--min-seq-id", str(identity),
        "--cov-mode", "0",
        "-c", "0.8",
        "-s", "7"
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Failed to cluster proteins: {e.stderr}\n")
        sys.exit(1)
    
    rep_seq_file = f"{output_prefix}_rep_seq.fasta"
    print(f"Representative sequences saved to {rep_seq_file}")
    return rep_seq_file


def run_mmseqs2_search(protein_db: str, genome_file: str, output_dir: str, evalue: float = 1e-5) -> str:
    """
    Run mmseqs2 search (tblastn-like) against genome.
    
    Parameters
    ----------
    protein_db : str
        Path to clustered protein sequences
    genome_file : str
        Path to target genome FASTA file
    output_dir : str
        Directory for output files
    evalue : float
        E-value threshold
        
    Returns
    -------
    str
        Path to search results file
    """
    mmseqs_dir = os.path.join(output_dir, "mmseqs_dbs")
    results_dir = os.path.join(output_dir, "mmseqs_results")
    os.makedirs(mmseqs_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    # Create databases
    protein_db_path = os.path.join(mmseqs_dir, "protein_db")
    genome_db_path = os.path.join(mmseqs_dir, "genome_db")
    
    print("Creating mmseqs2 databases...")
    subprocess.run(["mmseqs", "createdb", protein_db, protein_db_path], check=True)
    subprocess.run(["mmseqs", "createdb", genome_file, genome_db_path], check=True)
    
    # Run search
    search_result = os.path.join(results_dir, "search_result")
    tmp_dir = os.path.join(mmseqs_dir, "tmp")
    
    print("Running mmseqs2 search...")
    search_cmd = [
        "mmseqs", "search",
        protein_db_path, genome_db_path, search_result, tmp_dir,
        "--translation-table", "11",
        "-e", str(evalue)
    ]
    subprocess.run(search_cmd, check=True)
    
    # Convert to readable format
    output_file = os.path.join(results_dir, "search_results.txt")
    convert_cmd = [
        "mmseqs", "convertalis",
        protein_db_path, genome_db_path, search_result, output_file,
        "--format-output", "query,qstart,qend,qcov,target,tstart,tend,tframe,evalue,pident,alnlen,qlen"
    ]
    subprocess.run(convert_cmd, check=True)
    
    print(f"Search results saved to {output_file}")
    return output_file


def run_pseudomancer_pipeline(genus: str, genome_file: str, output_dir: str, evalue: float = 1e-5):
    """
    Run the complete pseudomancer pipeline with mmseqs2.
    
    Parameters
    ----------
    genus : str
        Genus name for protein download
    genome_file : str
        Path to target genome file
    output_dir : str
        Output directory
    evalue : float
        E-value threshold for searches
    """
    check_dependencies()
    
    # Download and merge proteins
    protein_file = download_ncbi_assemblies(genus, output_dir)
    
    # Cluster proteins
    clustered_proteins = cluster_proteins_mmseqs2(protein_file, output_dir)
    
    # Run search
    results_file = run_mmseqs2_search(clustered_proteins, genome_file, output_dir, evalue)
    
    print(f"Pipeline completed! Results in {results_file}")
    return results_file
