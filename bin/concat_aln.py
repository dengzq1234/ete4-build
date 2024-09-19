#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
from glob import glob
import sys
from argparse import ArgumentParser

def load_target_species(taxa_file):
    """Load expected taxa from a file."""
    target_species = set()
    with open(taxa_file, 'r') as f:
        for line in f:
            target_species.add(line.strip())
    return target_species

def detect_sequence_type(seq):
    """Detect whether a sequence is DNA or protein."""
    dna_chars = set("ACGTURYKMSWBDHVN")
    aa_chars = set("ACDEFGHIKLMNPQRSTVWYBXZ")
    
    seq = seq.upper().replace("-", "")  # Remove gaps and convert to uppercase
    
    # Count number of nucleotide and amino acid characters
    dna_count = sum(1 for char in seq if char in dna_chars)
    aa_count = sum(1 for char in seq if char in aa_chars)
    
    # Return 'DNA' or 'PROTEIN' based on which count is higher
    return "DNA" if dna_count > aa_count else "PROTEIN"

def process_alignment_files(alignment_files, target_species, separator, field=0):
    """Process the alignment files and create a concatenated supermatrix and partition data."""
    concat = defaultdict(str)
    partitions = []
    alignment_start = 1

    # Process each alignment file
    for algfile in alignment_files:
        visited_species = set()
        seq_lengths = []
        sequence_type = None  # To store whether it's DNA or AA

        # Read the alignment file
        try:
            with open(algfile, 'r') as f:
                for record in SeqIO.parse(f, format='fasta'):
                    spcode = record.id.split(separator)[field]
                    visited_species.add(spcode)
                    seq_lengths.append(len(record.seq))
                    concat[spcode] += str(record.seq)

                    # Detect sequence type if not already detected
                    if not sequence_type:
                        sequence_type = detect_sequence_type(str(record.seq))
        except FileNotFoundError:
            sys.exit(f"Error: Alignment file {algfile} not found.")
        except Exception as e:
            sys.exit(f"Error reading {algfile}: {str(e)}")

        # Ensure all sequences in the alignment have the same length
        if len(set(seq_lengths)) != 1:
            sys.exit(f"Error: Sequences in {algfile} have varying lengths.")

        # Calculate partition end
        alignment_end = alignment_start + seq_lengths[0] - 1
        partitions.append(f"{sequence_type}, gene_{algfile} = {alignment_start}-{alignment_end}")
        alignment_start = alignment_end + 1

        # Add gaps for missing species
        gap_line = '-' * seq_lengths[0]
        for missing_sp in target_species - visited_species:
            concat[missing_sp] += gap_line

    return concat, partitions

def write_supermatrix(output_file, concat):
    """Write the concatenated supermatrix to the output file."""
    try:
        with open(output_file, 'w') as output:
            for sp, seq in concat.items():
                output.write(f">{sp}\n{seq}\n")
    except Exception as e:
        sys.exit(f"Error writing to output file: {str(e)}")

def write_partition_file(partition_file, partitions):
    """Write the partition file if the flag is set."""
    try:
        with open(partition_file, 'w') as partition_out:
            for partition in partitions:
                partition_out.write(partition + '\n')
    except Exception as e:
        sys.exit(f"Error writing to partition file: {str(e)}")

def main():
    # Argument parser setup
    parser = ArgumentParser(description='Concatenate multiple alignments into a supermatrix and generate a partition file.')

    parser.add_argument('-a', dest='alignment_files', help='A list of alignment files in FASTA format (supports wildcards).', nargs="+", required=True)
    parser.add_argument('--spname-delimiter', dest='separator', help='Separator for species code in sequence names.', default='|')
    parser.add_argument('--spname-field', dest='field', help='Field number for species code in sequence names.', default='1')
    parser.add_argument('--taxa', dest='target_taxa', help='A file containing a list of expected taxa.', required=True)
    parser.add_argument('-o', dest='output_file', help='Output file to store the supermatrix.', required=True)
    parser.add_argument('-p', dest='partition_file', help='(Optional) Output file to store the partition information.', required=False)

    args = parser.parse_args()

    # Expand wildcards in alignment files using glob
    expanded_alignment_files = []
    for pattern in args.alignment_files:
        expanded_alignment_files.extend(glob(pattern))

    # Ensure at least one file is found
    if not expanded_alignment_files:
        sys.exit("Error: No alignment files found.")

    # Load expected taxa
    target_species = load_target_species(args.target_taxa)

    # Process the alignment files
    concat, partitions = process_alignment_files(expanded_alignment_files, target_species, args.separator, int(args.field))

    # Output the concatenated supermatrix
    write_supermatrix(args.output_file, concat)
    
    # If partition file is specified, write it
    if args.partition_file:
        write_partition_file(args.partition_file, partitions)
        print(f"Partition file written to {args.partition_file}")
    
    print(f"Supermatrix alignment written to {args.output_file}")

if __name__ == "__main__":
    main()
