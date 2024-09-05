#!/usr/bin/python3
import os
import sys
import argparse
from Bio import SeqIO

# Argument parser setup
parser = argparse.ArgumentParser(description="Convert between FASTA and PHYLIP formats")
subparsers = parser.add_subparsers(dest="command", help="Sub-command to specify conversion direction")

# Subcommand: fasta2phylip
fasta2phylip_parser = subparsers.add_parser('fasta2phylip', help="Convert FASTA to PHYLIP")
fasta2phylip_parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
fasta2phylip_parser.add_argument("-o", "--output", required=True, help="Output PHYLIP file")

# Subcommand: phylip2fasta
phylip2fasta_parser = subparsers.add_parser('phylip2fasta', help="Convert PHYLIP to FASTA")
phylip2fasta_parser.add_argument("-i", "--input", required=True, help="Input PHYLIP file")
phylip2fasta_parser.add_argument("-o", "--output", required=True, help="Output FASTA file")

args = parser.parse_args()

# Determine conversion direction and set input/output formats accordingly
if args.command == "fasta2phylip":
    input_format = "fasta"
    output_format = "phylip"
elif args.command == "phylip2fasta":
    input_format = "phylip"
    output_format = "fasta"
else:
    parser.print_help()
    sys.exit(1)

# Input and output file paths
input_file = args.input
output_file = args.output

# Convert between formats
records = SeqIO.parse(input_file, input_format)
count = SeqIO.write(records, output_file, output_format)