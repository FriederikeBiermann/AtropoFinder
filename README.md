# AtropoFinder
# Corefiner

## Overview
Corefiner makes identifying RiPP cores easy. It reads a specified FASTA file, fetches related genomic data from NCBI, and identifies regions of interest around bait enzymes to find precursor genes and core peptides.

## Prerequisites
Before you start, ensure you have the following software installed:
- Python 3.x
- Biopython
- pandas

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/corefiner.git
   cd corefiner
2. **Install dependencies:**
   ```bash
   pip install biopython pandas
## Usage
### Command-line Arguments
    -i, --input (required): Specifies the input FASTA file.
    -b, --boundary: The window around the bait enzyme in bp in which to look for precursor genes. Default is 3000 bp.
    -d, --dynamic_core_detection: Enables dynamic core detection (default: True).
    -e, --email (required): Email account to use for NCBI Entrez.
    -o, --output: Specifies the output directory. Default is Output/.
### Example code
  ```bash
  python corefiner.py -i input.fasta -b 3000 -d True -e your_email@example.com -o Output/
```
### Output Files
  info_precursor_peptides.csv: Information on the identified precursor peptides.
  Genbank_files/: Genbank files of identified regions.
  open_reading_frames.fasta: FASTA file of open reading frames.
  core_peptides.fasta: FASTA file of core peptides.
