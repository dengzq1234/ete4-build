#!/usr/bin/python3
import os
import sys
import argparse
from Bio import SeqIO

# Argument parser setup
parser = argparse.ArgumentParser(description="Convert between FASTA, PHYLIP, and NEXUS formats")
subparsers = parser.add_subparsers(dest="command", help="Sub-command to specify conversion direction")

# Subcommand: fasta2phylip
fasta2phylip_parser = subparsers.add_parser('fasta2phylip', help="Convert FASTA to PHYLIP")
fasta2phylip_parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
fasta2phylip_parser.add_argument("-o", "--output", required=True, help="Output PHYLIP file")

# Subcommand: phylip2fasta
phylip2fasta_parser = subparsers.add_parser('phylip2fasta', help="Convert PHYLIP to FASTA")
phylip2fasta_parser.add_argument("-i", "--input", required=True, help="Input PHYLIP file")
phylip2fasta_parser.add_argument("-o", "--output", required=True, help="Output FASTA file")

# Subcommand: fasta2nexus
fasta2nexus_parser = subparsers.add_parser('fasta2nexus', help="Convert FASTA to NEXUS")
fasta2nexus_parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
fasta2nexus_parser.add_argument("-o", "--output", required=True, help="Output NEXUS file")

args = parser.parse_args()

def detect_sequence_type(sequence):
    """
    Detect the type of the sequence as DNA, RNA, or protein based on the characters present.
    """
    dna_bases = set("-ACGTN")  # Standard DNA bases with ambiguous N
    rna_bases = set("-ACGUN")  # RNA bases with U instead of T
    protein_bases = set("-ACDEFGHIKLMNPQRSTVWYBXZ")  # Common protein bases including ambiguous X, B, Z
    seq_set = set(sequence.upper())

    # Check if the sequence fits DNA, RNA, or protein criteria
    if seq_set.issubset(dna_bases):
        return "DNA"
    elif seq_set.issubset(rna_bases):
        return "RNA"
    elif seq_set.issubset(protein_bases):
        return "protein"
    else:
        return

# Determine conversion direction and set input/output formats accordingly
if args.command == "fasta2phylip":
    input_format = "fasta"
    output_format = "phylip"
elif args.command == "phylip2fasta":
    input_format = "phylip"
    output_format = "fasta"
elif args.command == "fasta2nexus":
    input_format = "fasta"
    output_format = "nexus"
else:
    parser.print_help()
    sys.exit(1)

# Input and output file paths
input_file = args.input
output_file = args.output

# Convert between formats
records = SeqIO.parse(input_file, input_format)
# Convert between formats with appropriate settings for each conversion
if args.command == "fasta2nexus":
    record = next(records)
    molecule_type = detect_sequence_type(record.seq)
    records = SeqIO.convert(input_file, input_format, output_file, output_format, molecule_type)
    count = SeqIO.convert(input_file, input_format, output_file, output_format, molecule_type)
else:
    # Convert for fasta2phylip and phylip2fasta
    records = SeqIO.parse(input_file, input_format)
    count = SeqIO.write(records, output_file, output_format)

print(f"Converted {count} records from {input_file} to {output_file} in {output_format} format.")
# except Exception as e:
#     print(f"Error during conversion: {e}")
#     sys.exit(1)
