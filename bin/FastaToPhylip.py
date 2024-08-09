#!/usr/bin/python3
import os
import sys
import argparse
from Bio import SeqIO

# Argument parser setup
parser = argparse.ArgumentParser(description="Convert FASTA files to PHYLIP format")
parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
parser.add_argument("-o", "--output", required=True, help="Output PHYLIP file")

args = parser.parse_args()

# Input and output file paths
input_file = args.input
output_file = args.output

# Convert FASTA to PHYLIP
records = SeqIO.parse(input_file, "fasta")
count = SeqIO.write(records, output_file, "phylip")

print(f"Converted {count} records from {input_file} to {output_file}")