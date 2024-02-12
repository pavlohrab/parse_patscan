# Usage
This script for parsing Patscan files in a format:

```
NC_003888.3:[166177,166228]  : GCGGCACA CGCC TCGGCACA CGGC GCCGCACA ACGACGCGCGGA GCGGCACT
NC_003888.3:[554490,554545]  : CCCGCACG CCGGAGTC CCCGCACG CCGGCGTC CCCGCGCG CCCTCCTC CCCGCGCG
NC_003888.3:[794804,794850]  : CCCGCACC T CCCGCATC ACCCTCCTCG CCGGCACC TGCG CCGGCACC
NC_003888.3:[1260749,1260798]: CCCGCTCA GACA CCCGCTCA GACATCCCGA CCCGCTCA GACA CCCGCTCA
NC_003888.3:[1829090,1829124]: CCCGCACC G CCCGCACC G CCCGCACC G CCCGCACC
NC_003888.3:[1831031,1831080]: GCCGCACC GGGGCCGCCC GCCGCGCC GCGA GCCGCACC CGCA CCCGCACC
NC_003888.3:[3683264,3683313]: CCCGCACC GACGCCGCCCGTG CCCTCACC GGCG CCCGCGCC C TCCGCACC
```

Given patscan file and Genbank file, the script will retrieve CDS sequences, based in genome coordinates from patscan file. If the subsequence from Patscan is fully within range of CDS, if will be classified accordingly and information about CDS will be added to the csv file. Otherwise it will be classified as intergenic 

Optionally user can write faa and fna files with unique sequences.

```{bash}
usage: parse_patscan.py [-h] [-g GENBANK_FILE] [-p PATSCAN_FILE] [-o OUTPUT_FILE] [--write_faa] [--write_fna] [-v]

This script takes a genbank file and a patscan file and returns a csv file with the location of the patscan sequences in the genbank file

options:
  -h, --help            show this help message and exit
  -g GENBANK_FILE, --genbank_file GENBANK_FILE
                        The genbank file
  -p PATSCAN_FILE, --patscan_file PATSCAN_FILE
                        The patscan file
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        The output file
  --write_faa           Write a fasta file with the protein sequences
  --write_fna           Write a fasta file with the nucleotide sequences
  -v, --version         show program's version number and exit

Version 0.1.0 - (C) Pavlo Hrab
```

# Installation

You need biopython and pandas in an environment. You ca install them with conda/mamba/micromamba:

```{bash}
conda create -n patscan biopython pandas
conda activate patscan
```
