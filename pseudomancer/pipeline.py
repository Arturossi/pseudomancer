"""
Main pipeline orchestration for pseudomancer.
"""

from .dependencies import check_dependencies
from .download import download_ncbi_assemblies
from .clustering import cluster_proteins_mmseqs2
from .orfs import identify_orfs
from .search import run_mmseqs2_search
from .candidates import process_mmseqs2_results
from .overlap import overlap_with_getorf


def run_pseudomancer_pipeline(taxon: str, genome_file: str, output_dir: str, evalue: float = 1e-5):
    """
    Run the complete pseudomancer pipeline with mmseqs2.
    
    Parameters
    ----------
    taxon : str
        Taxon name for protein download
    genome_file : str
        Path to target genome file
    output_dir : str
        Output directory
    evalue : float
        E-value threshold for searches
    """
    check_dependencies()

    # ORF identification in reference genome
    orf_file, orf_gff_file = identify_orfs(genome_file, output_dir)

    # Download and merge proteins
    protein_file = download_ncbi_assemblies(taxon, output_dir)

    # Cluster proteins
    clustered_proteins = cluster_proteins_mmseqs2(protein_file, output_dir)

    # Run search
    results_file = run_mmseqs2_search(clustered_proteins, genome_file, output_dir, evalue)

    # Process mmseqs2 results to generate pseudogene candidates GFF
    pseudogene_gff = process_mmseqs2_results(results_file, output_dir, genome_file)

    # Assign clustered tblastn hits to ORFs
    mappings_file = overlap_with_getorf(pseudogene_gff, orf_gff_file, output_dir)
    
    print(f"Pipeline completed! ORF mappings in {mappings_file}")
    return mappings_file