#!/usr/bin/python3
from Bio import SeqIO
import os
import sys

fasta_files = sys.argv[1:]
for file in fasta_files:
    fasta_name = os.path.basename(file).replace(".clean.alg", "")
    records = SeqIO.parse(file, "fasta")
    phylip_seq = f'{fasta_name}.clean.alg.phylip'
    count = SeqIO.write(records, phylip_seq, "phylip")
    
    # # Print a separator if there are multiple files
    # if len(fasta_files) > 1:
    #     print("\n" + "="*80 + "\n")