# AtropoFinder
## Overview

This Jupyter Notebook provides an analysis pipeline for identifying and classifying tryptorubin-like peptides using sequence motifs.
## Prerequisites

    Python 3
    Biopython 1.76
    Pandas
    Numpy
    Scikit-learn 1.0.2
    Imbalanced-learn 0.9.0

## Installation

To install the required packages, you can use pip:
   ```bash
   pip install biopython==1.76 pandas numpy scikit-learn==1.0.2 imbalanced-learn==0.9.0

   ```
## Usage

    Prepare Input Files:
        Place your FASTA files in the Fasta directory.
        Ensure the files are correctly named and paths are set in the notebook TerP450-TLP Version.ipynb.

    Set Parameters:
        Adjust the LENGTH variable to set the length of the sequence motif you want to test.
        Update foldernameoutput to the desired output directory.

    Run the Notebook:
        Execute the cells in the notebook TerP450-TLP Version.ipynb to process the train the classifier.
        To run the pretrained classifier on you own sequences, align them using Clustal W version 1.2.3 with the default settings to the reference protein cytochrome P450 with known annotation of functional regions from Mycobacterium tuberculosis H2102 (GenBank: KBE51585.1) in a multiple sequence alignment. Fragment the alignment at positions 92, 192, 275, and 395. Add the respecting files in the cell " Analyse sequences with fragment fastas" and run the cell.
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
