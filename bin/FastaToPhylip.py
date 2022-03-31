#!/usr/bin/python3
from Bio import SeqIO
import sys

fasta_seq = sys.argv[1]
phylip_seq = 'clean.alg.phylip'

records = SeqIO.parse(fasta_seq, "fasta")
count = SeqIO.write(records, phylip_seq, "phylip")