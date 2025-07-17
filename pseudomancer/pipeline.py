import subprocess
import sys
import os
import shutil
import tempfile
import requests
import pandas as pd
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any
import zipfile
import re
from collections import defaultdict
import math


def check_dependencies():
    """Check if required dependencies (mmseqs2, datasets, getorf) are available."""
    dependencies = ["mmseqs", "datasets", "getorf"]
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
            "- getorf (EMBOSS): conda install -c bioconda emboss\n"
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
        "datasets",
        "download",
        "genome",
        "taxon",
        genus,
        "--annotated",
        "--reference",
        "--assembly-level",
        "complete",
        "--include",
        "genome,protein,gff3,gtf,seq-report",
        "--filename",
        os.path.join(assemblies_dir, f"{genus.lower()}_dataset.zip"),
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Failed to download assemblies: {e.stderr}\n")
        sys.exit(1)

    # Extract and merge protein files
    protein_files = []
    zip_path = os.path.join(assemblies_dir, f"{genus.lower()}_dataset.zip")

    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(assemblies_dir)

    # Find all protein.faa files
    for root, dirs, files in os.walk(assemblies_dir):
        for file in files:
            if file == "protein.faa":
                protein_files.append(os.path.join(root, file))

    # Merge protein files
    merged_proteins = os.path.join(output_dir, f"{genus.lower()}_proteins.faa")
    with open(merged_proteins, "w") as outfile:
        for protein_file in protein_files:
            with open(protein_file, "r") as infile:
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
        "mmseqs",
        "easy-cluster",
        protein_file,
        output_prefix,
        output_dir,
        "--min-seq-id",
        str(identity),
        "--cov-mode",
        "0",
        "-c",
        "0.8",
        "-s",
        "7",
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Failed to cluster proteins: {e.stderr}\n")
        sys.exit(1)

    rep_seq_file = f"{output_prefix}_rep_seq.fasta"
    print(f"Representative sequences saved to {rep_seq_file}")
    return rep_seq_file


def parse_getorf_header(header: str) -> Optional[dict]:
    """
    Parse getorf FASTA header to extract coordinates.

    Example headers:
    Forward: >NC_002677.1_1 [93 - 272] Mycobacterium leprae TN, complete sequence
    Reverse: >NC_002677.1_53466 [1611 - 1420] (REVERSE SENSE) Mycobacterium leprae TN, complete sequence
    """
    match = re.match(r"^(.+)_(\d+)", header)
    if not match:
        return None

    seq_id = match.group(1)
    orf_num = match.group(2)

    coord_match = re.search(r"\[(\d+) - (\d+)\]", header)
    if not coord_match:
        return None

    start = int(coord_match.group(1))
    end = int(coord_match.group(2))

    if "(REVERSE SENSE)" in header:
        strand = "-"
    else:
        strand = "+"

    if strand == "+":
        frame = ((start - 1) % 3) + 1
    else:
        frame = -((end - 1) % 3 + 1)

    return {
        "seq_id": seq_id,
        "orf_num": orf_num,
        "start": start,
        "end": end,
        "strand": strand,
        "frame": frame,
        "length": abs(end - start) + 1,
    }


def convert_fasta_to_gff(fasta_file: str, gff_file: str) -> None:
    """Convert getorf FASTA output to GFF format."""
    print(f"Converting {fasta_file} to GFF format...")

    gff_lines = ["##gff-version 3"]
    orf_count = 0

    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Parse header
                header = line.strip()[1:]  # Remove '>'
                orf_info = parse_getorf_header(header)
                if not orf_info:
                    print(f"Warning: Could not parse header: {header}")
                    continue

                orf_count += 1

                # Create GFF attributes
                orf_id = f"getorf_orf_{orf_info['orf_num']}"
                attributes = f"ID={orf_id};Name={orf_id};length={orf_info['length']};frame={orf_info['frame']}"

                # Create GFF line
                gff_line = "\t".join(
                    [
                        orf_info["seq_id"],  # seqid
                        "getorf",  # source
                        "ORF",  # type
                        str(orf_info["start"]),  # start
                        str(orf_info["end"]),  # end
                        ".",  # score
                        orf_info["strand"],  # strand
                        ".",  # phase
                        attributes,  # attributes
                    ]
                )

                gff_lines.append(gff_line)

    with open(gff_file, "w") as f:
        f.write("\n".join(gff_lines) + "\n")

    print(f"Converted {orf_count} ORFs to GFF format")
    print(f"GFF output written to {gff_file}")


def identify_orfs(genome_file: str, output_dir: str, minsize: int = 90, table: int = 11) -> Tuple[str, str]:
    """
    Identify open reading frames (ORFs) in the genome using getorf.

    Parameters
    ----------
    genome_file : str
        Path to genome FASTA file
    output_dir : str
        Directory for output files
    minsize : int
        Minimum ORF size in codons (default: 90)
    table : int
        Genetic code table (default: 11 for bacteria)

    Returns
    -------
    Tuple[str, str]
        Paths to ORF sequences file and GFF file
    """
    base_name = os.path.splitext(os.path.basename(genome_file))[0]
    orf_file = os.path.join(output_dir, f"{base_name}_orfs.fasta")
    gff_file = os.path.join(output_dir, f"{base_name}_orfs.gff")

    print(f"Running getorf ORF identification on {genome_file}...")

    cmd = [
        "getorf",
        "-sequence",
        genome_file,
        "-outseq",
        orf_file,
        "-find",
        "3",
        "-minsize",
        str(minsize),
        "-table",
        str(table),
        "-reverse",
        "yes",
        "-flanking",
        "0",
        "-methionine",
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Failed to identify ORFs: {e.stderr}\n")
        sys.exit(1)

    print(f"ORF sequences saved to {orf_file}")

    convert_fasta_to_gff(orf_file, gff_file)

    return orf_file, gff_file


def calculate_distance(start1: int, end1: int, start2: int, end2: int) -> int:
    """Calculate distance between two genomic regions."""
    return max(0, max(min(start1, end1), min(start2, end2)) - min(max(start1, end1), max(start2, end2)))


def are_adjacent(start1: int, end1: int, start2: int, end2: int, max_distance: int = 100) -> bool:
    """Check if two genomic regions are adjacent within max_distance."""
    return calculate_distance(start1, end1, start2, end2) <= max_distance


def group_adjacent_hits(hits_df: pd.DataFrame, max_distance: int = 100) -> pd.DataFrame:
    """Group adjacent hits within max_distance."""
    if len(hits_df) == 0:
        return hits_df
    
    # Sort by genomic position
    hits_df = hits_df.sort_values(by=['tstart', 'tend'])
    hits_df = hits_df.reset_index(drop=True)
    
    # Initialize clustering
    clusters = []
    hits_df['cluster'] = 0
    current_cluster = 1
    
    hits_df.loc[0, 'cluster'] = current_cluster
    
    for i in range(1, len(hits_df)):
        assigned = False
        
        # Check against all existing clusters
        for cluster_id in range(1, current_cluster + 1):
            cluster_hits = hits_df[hits_df['cluster'] == cluster_id]
            
            # Check if current hit is adjacent to any hit in this cluster
            for j in cluster_hits.index:
                if are_adjacent(hits_df.loc[i, 'tstart'], hits_df.loc[i, 'tend'], 
                              hits_df.loc[j, 'tstart'], hits_df.loc[j, 'tend'], max_distance):
                    hits_df.loc[i, 'cluster'] = cluster_id
                    assigned = True
                    break
            if assigned:
                break
        
        # If not assigned to any existing cluster, create new cluster
        if not assigned:
            current_cluster += 1
            hits_df.loc[i, 'cluster'] = current_cluster
    
    return hits_df


def create_cluster_summary(cluster_hits: pd.DataFrame, query_name: str, event_type: str, strand_name: str) -> Dict[str, Any]:
    """Create summary for a cluster of hits."""
    # Sort hits by genomic position
    cluster_hits = cluster_hits.sort_values(by=['tstart', 'tend'])
    
    # Calculate cluster boundaries
    cluster_start = min(cluster_hits['tstart'].min(), cluster_hits['tend'].min())
    cluster_end = max(cluster_hits['tstart'].max(), cluster_hits['tend'].max())
    
    # Calculate query boundaries
    query_start = cluster_hits['qstart'].min()
    query_end = cluster_hits['qend'].max()
    
    # Get reading frames involved
    frames = ','.join(map(str, sorted(cluster_hits['tframe'].unique())))
    
    # Calculate best e-value and mean percent identity
    best_evalue = cluster_hits['evalue'].min()
    mean_pident = cluster_hits['pident'].mean()
    
    # Create hit details string
    hit_details = []
    for _, hit in cluster_hits.iterrows():
        hit_detail = f"q{hit['qstart']}-{hit['qend']}:t{hit['tstart']}-{hit['tend']}:f{hit['tframe']}:e{hit['evalue']:.2e}"
        hit_details.append(hit_detail)
    
    return {
        'query': query_name,
        'target': cluster_hits['target'].iloc[0],
        'cluster_id': cluster_hits['cluster'].iloc[0],
        'event_type': event_type,
        'strand': strand_name,
        'num_hits': len(cluster_hits),
        'frames': frames,
        'query_start': query_start,
        'query_end': query_end,
        'query_span': query_end - query_start + 1,
        'target_start': cluster_start,
        'target_end': cluster_end,
        'target_span': cluster_end - cluster_start + 1,
        'best_evalue': best_evalue,
        'mean_pident': mean_pident,
        'hit_details': ';'.join(hit_details)
    }


def process_query(query_data: pd.DataFrame) -> List[Dict[str, Any]]:
    """Process a single query to identify frameshift and nonsense candidates."""
    query_name = query_data['query'].iloc[0]
    
    # If only one hit, return empty list
    if len(query_data) == 1:
        return []
    
    # Separate hits by strand
    forward_hits = query_data[query_data['tframe'] > 0]
    reverse_hits = query_data[query_data['tframe'] < 0]
    
    candidates = []
    
    # Process forward strand hits
    if len(forward_hits) > 1:
        forward_clustered = group_adjacent_hits(forward_hits)
        for cluster_id in forward_clustered['cluster'].unique():
            cluster_hits = forward_clustered[forward_clustered['cluster'] == cluster_id]
            
            if len(cluster_hits) == 1:
                continue  # Skip single-hit clusters
            
            # Check if cluster contains multiple reading frames
            unique_frames = cluster_hits['tframe'].abs().unique()
            
            if len(unique_frames) > 1:
                # Multiple reading frames - frameshift candidate
                candidate = create_cluster_summary(cluster_hits, query_name, 'frameshift', 'forward')
                candidates.append(candidate)
            else:
                # Same reading frame - nonsense candidate
                candidate = create_cluster_summary(cluster_hits, query_name, 'nonsense', 'forward')
                candidates.append(candidate)
    
    # Process reverse strand hits
    if len(reverse_hits) > 1:
        reverse_clustered = group_adjacent_hits(reverse_hits)
        for cluster_id in reverse_clustered['cluster'].unique():
            cluster_hits = reverse_clustered[reverse_clustered['cluster'] == cluster_id]
            
            if len(cluster_hits) == 1:
                continue  # Skip single-hit clusters
            
            # Check if cluster contains multiple reading frames
            unique_frames = cluster_hits['tframe'].abs().unique()
            
            if len(unique_frames) > 1:
                # Multiple reading frames - frameshift candidate
                candidate = create_cluster_summary(cluster_hits, query_name, 'frameshift', 'reverse')
                candidates.append(candidate)
            else:
                # Same reading frame - nonsense candidate
                candidate = create_cluster_summary(cluster_hits, query_name, 'nonsense', 'reverse')
                candidates.append(candidate)
    
    return candidates


def cluster_overlapping_candidates(candidates_df: pd.DataFrame, min_overlap: int = 1) -> pd.DataFrame:
    """Cluster overlapping genomic regions from different query proteins."""
    if len(candidates_df) == 0:
        return candidates_df
    
    # Determine strand from frames column
    candidates_df['strand_direction'] = candidates_df['frames'].apply(
        lambda x: 'reverse' if any(int(f) < 0 for f in str(x).split(',')) else 'forward'
    )
    
    # Group by target and strand
    grouped = candidates_df.groupby(['target', 'strand_direction'])
    
    merged_results = []
    
    for (target, strand_direction), group_df in grouped:
        if len(group_df) <= 1:
            # Single candidate - create merged entry
            for _, row in group_df.iterrows():
                merged_results.append({
                    'merged_id': f"merged_{target}_{strand_direction}_1",
                    'target': target,
                    'strand': strand_direction,
                    'target_start': row['target_start'],
                    'target_end': row['target_end'],
                    'target_span': row['target_span'],
                    'num_contributing_queries': 1,
                    'contributing_queries': row['query'],
                    'all_frames': row['frames'],
                    'best_evalue': row['best_evalue'],
                    'mean_pident': row['mean_pident'],
                    'max_pident': row['mean_pident'],
                    'total_hits': row['num_hits'],
                    'min_query_start': row['query_start'],
                    'max_query_end': row['query_end'],
                    'total_query_span': row['query_span'],
                    'event_type': row['event_type'],
                    'combined_hit_details': row['hit_details']
                })
            continue
        
        # Find overlapping candidates using simple approach
        group_df = group_df.reset_index(drop=True)
        visited = [False] * len(group_df)
        overlap_group = 0
        
        for i in range(len(group_df)):
            if visited[i]:
                continue
                
            overlap_group += 1
            current_group = [i]
            visited[i] = True
            
            # Find all overlapping candidates
            for j in range(i + 1, len(group_df)):
                if visited[j]:
                    continue
                    
                # Check if j overlaps with any in current_group
                for k in current_group:
                    overlap_start = max(group_df.iloc[j]['target_start'], group_df.iloc[k]['target_start'])
                    overlap_end = min(group_df.iloc[j]['target_end'], group_df.iloc[k]['target_end'])
                    
                    if overlap_start <= overlap_end and (overlap_end - overlap_start + 1) >= min_overlap:
                        current_group.append(j)
                        visited[j] = True
                        break
            
            # Create merged candidate for this overlap group
            group_candidates = group_df.iloc[current_group]
            
            merged_results.append({
                'merged_id': f"merged_{target}_{strand_direction}_{overlap_group}",
                'target': target,
                'strand': strand_direction,
                'target_start': group_candidates['target_start'].min(),
                'target_end': group_candidates['target_end'].max(),
                'target_span': group_candidates['target_end'].max() - group_candidates['target_start'].min() + 1,
                'num_contributing_queries': len(group_candidates),
                'contributing_queries': ';'.join(group_candidates['query'].unique()),
                'all_frames': ','.join(sorted(set(frame for frames in group_candidates['frames'] for frame in str(frames).split(',')))),
                'best_evalue': group_candidates['best_evalue'].min(),
                'mean_pident': group_candidates['mean_pident'].mean(),
                'max_pident': group_candidates['mean_pident'].max(),
                'total_hits': group_candidates['num_hits'].sum(),
                'min_query_start': group_candidates['query_start'].min(),
                'max_query_end': group_candidates['query_end'].max(),
                'total_query_span': group_candidates['query_span'].sum(),
                'event_type': ','.join(group_candidates['event_type'].unique()),
                'combined_hit_details': ';'.join(group_candidates['hit_details'])
            })
    
    return pd.DataFrame(merged_results)


def get_genome_info(genome_file: str) -> Tuple[str, int]:
    """Extract sequence ID and length from genome FASTA file."""
    with open(genome_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract sequence ID (first word after >)
                seq_id = line.strip()[1:].split()[0]
                break
        else:
            raise ValueError(f"No sequences found in {genome_file}")
    
    # Calculate total length
    total_length = 0
    with open(genome_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                total_length += len(line.strip())
    
    return seq_id, total_length


def export_candidates_to_gff(frameshift_df: pd.DataFrame, output_file: str, genome_file: str) -> None:
    """Export frameshift candidates to GFF3 format."""
    # Get genome information
    seq_id, genome_length = get_genome_info(genome_file)
    
    # Create GFF3 header
    gff_lines = [
        "##gff-version 3",
        f"##sequence-region {seq_id} 1 {genome_length}"
    ]
    
    # Process frameshift candidates
    for i, row in frameshift_df.iterrows():
        # Create display name
        gene_display = row['contributing_queries'].split(';')[0].replace('.', '_')
        
        # Create attributes string
        attributes = [
            f"ID={row['merged_id']}",
            f"Name={gene_display}",
            f"Note=Merged frameshift candidate with {row['num_contributing_queries']} contributing queries, {row['total_hits']} total hits in frames {row['all_frames']}",
            f"merged_id={row['merged_id']}",
            f"num_contributing_queries={row['num_contributing_queries']}",
            f"contributing_queries={row['contributing_queries']}",
            f"total_hits={row['total_hits']}",
            f"all_frames={row['all_frames']}",
            f"total_query_span={row['total_query_span']}",
            f"best_evalue={row['best_evalue']:.2e}",
            f"mean_pident={row['mean_pident']:.1f}",
            f"max_pident={row['max_pident']:.1f}",
            f"event_type={row['event_type']}",
            f"hit_details={row['combined_hit_details']}"
        ]
        
        # Create GFF line
        gff_line = '\t'.join([
            row['target'],          # seqid
            'MMseqs2',             # source
            'pseudogene',          # type
            str(int(row['target_start'])),  # start
            str(int(row['target_end'])),    # end
            str(round(-math.log10(row['best_evalue']), 1)),  # score
            '+' if row['strand'] == 'forward' else '-',  # strand
            '.',                   # phase
            ';'.join(attributes)   # attributes
        ])
        
        gff_lines.append(gff_line)
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(gff_lines) + '\n')
    
    print(f"GFF3 file written to: {output_file}")
    print(f"Features exported: {len(frameshift_df)} merged frameshift candidates")


def process_mmseqs2_results(results_file: str, output_dir: str, genome_file: str) -> str:
    """Process mmseqs2 results to identify pseudogene candidates and generate GFF."""
    print(f"Processing mmseqs2 results from {results_file}")
    
    # Read mmseqs2 results
    column_names = ['query', 'qstart', 'qend', 'qcov', 'target', 'tstart', 'tend', 'tframe', 'evalue', 'pident', 'alnlen', 'qlen']
    mmseqs_results = pd.read_csv(results_file, sep='\t', names=column_names)
    
    print(f"Loaded {len(mmseqs_results)} mmseqs2 hits")
    
    # Process each query separately
    all_candidates = []
    
    for query in mmseqs_results['query'].unique():
        query_data = mmseqs_results[mmseqs_results['query'] == query]
        candidates = process_query(query_data)
        all_candidates.extend(candidates)
    
    print(f"Found {len(all_candidates)} individual candidates")
    
    # Convert to DataFrame
    candidates_df = pd.DataFrame(all_candidates)
    
    if len(candidates_df) == 0:
        print("No candidates found")
        return ""
    
    # Filter for frameshift candidates only
    frameshift_candidates = candidates_df[candidates_df['event_type'] == 'frameshift']
    
    print(f"Found {len(frameshift_candidates)} frameshift candidates")
    
    # Remove duplicates based on target_start
    frameshift_candidates = frameshift_candidates.drop_duplicates(subset=['target_start'])
    
    # Cluster overlapping candidates
    frameshift_merged = cluster_overlapping_candidates(frameshift_candidates)
    
    print(f"Merged into {len(frameshift_merged)} clustered candidates")
    
    # Generate GFF file
    gff_file = os.path.join(output_dir, "pseudogene_candidates_annotated_merged.gff")
    export_candidates_to_gff(frameshift_merged, gff_file, genome_file)
    
    return gff_file


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

    protein_db_path = os.path.join(mmseqs_dir, "protein_db")
    genome_db_path = os.path.join(mmseqs_dir, "genome_db")

    print("Creating mmseqs2 databases...")
    subprocess.run(["mmseqs", "createdb", protein_db, protein_db_path], check=True)
    subprocess.run(["mmseqs", "createdb", genome_file, genome_db_path], check=True)

    search_result = os.path.join(results_dir, "search_result")
    tmp_dir = os.path.join(mmseqs_dir, "tmp")

    print("Running mmseqs2 search...")
    search_cmd = [
        "mmseqs",
        "search",
        protein_db_path,
        genome_db_path,
        search_result,
        tmp_dir,
        "--translation-table",
        "11",
        "-e",
        str(evalue),
    ]
    subprocess.run(search_cmd, check=True)

    # Convert to readable format
    output_file = os.path.join(results_dir, "search_results.txt")
    convert_cmd = [
        "mmseqs",
        "convertalis",
        protein_db_path,
        genome_db_path,
        search_result,
        output_file,
        "--format-output",
        "query,qstart,qend,qcov,target,tstart,tend,tframe,evalue,pident,alnlen,qlen",
    ]
    subprocess.run(convert_cmd, check=True)

    print(f"Search results saved to {output_file}")
    return output_file


def parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GFF attributes string into a dictionary."""
    attributes = {}
    
    # Handle hit_details specially since it contains semicolons
    if "hit_details=" in attr_string:
        before_hit_details, rest = attr_string.split("hit_details=", 1)
        
        # Find where hit_details ends (look for next attribute with =)
        hit_details_end = -1
        for i, char in enumerate(rest):
            if char == ";" and "=" in rest[i+1:i+50]:
                hit_details_end = i
                break
        
        if hit_details_end == -1:
            hit_details_value = rest
            after_hit_details = ""
        else:
            hit_details_value = rest[:hit_details_end]
            after_hit_details = rest[hit_details_end+1:]
        
        # Parse before and after normally
        for part in [before_hit_details, after_hit_details]:
            if part:
                for item in part.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        attributes[key] = value
        
        attributes["hit_details"] = hit_details_value
    else:
        # Normal parsing if no hit_details
        for item in attr_string.split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                attributes[key] = value
    
    return attributes


def parse_hit_details(hit_details: str) -> List[Dict[str, Any]]:
    """Parse hit_details string to extract individual hit information."""
    hits = []
    for hit in hit_details.split(";"):
        if not hit.strip():
            continue
        
        # Parse: q146-189:t2913374-2913508:f1:e8.38e-09
        parts = hit.split(":")
        if len(parts) >= 4:
            try:
                q_coords = parts[0].replace("q", "")
                t_coords = parts[1].replace("t", "")
                mmseqs_frame = int(parts[2].replace("f", ""))
                evalue = parts[3]
                if evalue.startswith("e"):
                    evalue = evalue[1:]
                
                t_start, t_end = map(int, t_coords.split("-"))
                
                # Determine strand from frame
                strand = "+" if mmseqs_frame > 0 else "-"
                
                # Convert MMseqs2 frame to biological frame
                if mmseqs_frame > 0:
                    biological_frame = mmseqs_frame + 1
                    if biological_frame > 3:
                        biological_frame = 1
                else:
                    biological_frame = mmseqs_frame + 1
                    if biological_frame == 0:
                        biological_frame = -3
                
                hits.append({
                    "q_coords": q_coords,
                    "t_start": t_start,
                    "t_end": t_end,
                    "frame": biological_frame,
                    "strand": strand,
                    "evalue": evalue,
                })
            except (ValueError, IndexError):
                continue
    return hits


def find_overlapping_orfs(hit: Dict[str, Any], orfs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Find ORFs that overlap with a hit in the same strand and frame."""
    overlapping = []
    
    hit_start = min(hit["t_start"], hit["t_end"])
    hit_end = max(hit["t_start"], hit["t_end"])
    
    for orf in orfs:
        # Check strand and frame match
        if orf["strand"] != hit["strand"] or orf["frame"] != hit["frame"]:
            continue
        
        # Check for overlap
        orf_start = min(orf["start"], orf["end"])
        orf_end = max(orf["start"], orf["end"])
        
        overlap_start = max(hit_start, orf_start)
        overlap_end = min(hit_end, orf_end)
        
        if overlap_start <= overlap_end:
            overlap_length = overlap_end - overlap_start + 1
            overlapping.append({
                "orf": orf,
                "overlap_length": overlap_length,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
            })
    
    return overlapping


def parse_getorf_gff(getorf_gff: str) -> List[Dict[str, Any]]:
    """Parse getorf GFF output to extract ORF information."""
    orfs = []
    
    with open(getorf_gff, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            
            seqid = fields[0]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            # Only process ORF features
            if feature_type.lower() in ["orf", "cds"]:
                # Extract ORF ID from attributes
                orf_id = None
                for attr in attributes.split(";"):
                    if attr.startswith("ID="):
                        orf_id = attr.split("=", 1)[1]
                        break
                
                if not orf_id:
                    orf_id = f"orf_{seqid}_{start}_{end}_{strand}"
                
                # Calculate reading frame from coordinates and strand
                if strand == "+":
                    frame = ((start - 1) % 3) + 1
                else:
                    frame = -((end - 1) % 3 + 1)
                
                orfs.append({
                    "orf_id": orf_id,
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "frame": frame,
                    "length": end - start + 1,
                    "attributes": attributes,
                })
    
    return orfs


def overlap_with_getorf(pseudogene_gff: str, getorf_gff: str, output_dir: str) -> str:
    """Find overlapping ORFs for pseudogene candidates."""
    print(f"Finding ORF overlaps between {pseudogene_gff} and {getorf_gff}")
    
    # Parse getorf ORFs
    orfs = parse_getorf_gff(getorf_gff)
    print(f"Parsed {len(orfs)} ORFs from getorf output")
    
    # Parse pseudogene clusters
    clusters = []
    with open(pseudogene_gff, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            
            seqid = fields[0]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = parse_gff_attributes(fields[8])
            
            # Only process pseudogene features with hit_details
            if feature_type == "pseudogene" and "hit_details" in attributes:
                clusters.append({
                    "seqid": seqid,
                    "cluster_start": start,
                    "cluster_end": end,
                    "cluster_strand": strand,
                    "cluster_id": attributes.get("ID", f"cluster_{len(clusters)}"),
                    "cluster_name": attributes.get("Name", "Unknown"),
                    "hit_details": attributes["hit_details"],
                    "attributes": attributes,
                })
    
    print(f"Found {len(clusters)} pseudogene clusters")
    
    # Process each cluster
    results = []
    
    for cluster in clusters:
        # Parse hits from this cluster
        hits = parse_hit_details(cluster["hit_details"])
        
        # For each hit, find overlapping ORFs
        for hit_idx, hit in enumerate(hits):
            hit_id = f"{cluster['cluster_id']}_hit_{hit_idx + 1}"
            
            # Find overlapping ORFs in same strand and frame
            overlapping_orfs = find_overlapping_orfs(hit, orfs)
            
            if overlapping_orfs:
                # If multiple ORFs overlap, choose the longest one
                best_orf = max(overlapping_orfs, key=lambda x: x["orf"]["length"])
                
                results.append({
                    "cluster_id": cluster["cluster_id"],
                    "cluster_name": cluster["cluster_name"],
                    "hit_id": hit_id,
                    "hit_start": hit["t_start"],
                    "hit_end": hit["t_end"],
                    "hit_frame": hit["frame"],
                    "hit_strand": hit["strand"],
                    "hit_evalue": hit["evalue"],
                    "orf_id": best_orf["orf"]["orf_id"],
                    "orf_start": best_orf["orf"]["start"],
                    "orf_end": best_orf["orf"]["end"],
                    "orf_frame": best_orf["orf"]["frame"],
                    "orf_strand": best_orf["orf"]["strand"],
                    "orf_length": best_orf["orf"]["length"],
                    "overlap_length": best_orf["overlap_length"],
                    "overlap_start": best_orf["overlap_start"],
                    "overlap_end": best_orf["overlap_end"],
                    "num_overlapping_orfs": len(overlapping_orfs),
                })
            else:
                # Still record it but with no ORF info
                results.append({
                    "cluster_id": cluster["cluster_id"],
                    "cluster_name": cluster["cluster_name"],
                    "hit_id": hit_id,
                    "hit_start": hit["t_start"],
                    "hit_end": hit["t_end"],
                    "hit_frame": hit["frame"],
                    "hit_strand": hit["strand"],
                    "hit_evalue": hit["evalue"],
                    "orf_id": None,
                    "orf_start": None,
                    "orf_end": None,
                    "orf_frame": None,
                    "orf_strand": None,
                    "orf_length": None,
                    "overlap_length": None,
                    "overlap_start": None,
                    "overlap_end": None,
                    "num_overlapping_orfs": 0,
                })
    
    # Convert to DataFrame and save TSV
    df = pd.DataFrame(results)
    output_file = os.path.join(output_dir, "cluster_orf_mappings.tsv")
    df.to_csv(output_file, sep="\t", index=False)
    
    print(f"Results saved to {output_file}")
    print(f"Total hit-ORF mappings: {len(results)}")
    print(f"Hits with matching ORFs: {len(df[df['orf_id'].notna()])}")
    print(f"Hits without matching ORFs: {len(df[df['orf_id'].isna()])}")
    
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

    # ORF identification in reference genome
    orf_file, orf_gff_file = identify_orfs(genome_file, output_dir)

    # Download and merge proteins
    protein_file = download_ncbi_assemblies(genus, output_dir)

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
