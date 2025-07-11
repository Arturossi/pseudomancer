# pseudomancer

A command line tool for reconstructing pseudogenes in prokaryotic genomes using mmseqs2.

## Installation

### Prerequisites

Before using pseudomancer, you need to install these external dependencies:

1. **mmseqs2** - For protein clustering and homology searches
   ```bash
   # Install via conda (recommended)
   conda install -c conda-forge -c bioconda mmseqs2
   
   # Or see: https://github.com/soedinglab/MMseqs2
   ```

2. **NCBI datasets** - For downloading reference genomes
   ```bash
   # Download the appropriate binary for your system
   # See: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/
   
   # Example for Linux:
   wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
   chmod +x datasets
   # Move to a directory in your PATH
   ```

### Python Package

```bash
git clone <this-repo>
cd pseudomancer
python -m venv env
source env/bin/activate  # On Windows: env\Scripts\activate
pip install -r requirements.txt
```

## Usage

```bash
python -m pseudomancer --genus Mycobacterium --genome target_genome.fasta --out_dir results/
```

### Parameters

- `--genus`: Genus name for downloading reference proteins from NCBI RefSeq
- `--genome`: Target genome FASTA file to search for pseudogenes  
- `--out_dir`: Output directory for results
- `--evalue`: E-value threshold (default: 1e-5)

## How it works

1. Downloads all complete, annotated genomes for the specified genus from NCBI RefSeq
2. Extracts and merges protein sequences from all assemblies
3. Clusters proteins at 99% identity using mmseqs2 to create a non-redundant dataset
4. Searches the clustered proteins against your target genome using mmseqs2 (tblastn-like search)
5. Outputs results in tabular format with alignment statistics

## Note

The datasets download approach is experimental - for large genera, you may need to download assembly lists manually as shown in `test_data/new_pipeline.sh`.
