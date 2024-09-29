#!/usr/bin/env python3
import sys
import argparse
from ete4 import Tree

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Modify leaf names of phylogenetic trees based on specified delimiter and field.")
    parser.add_argument('files', nargs='+', help="Input tree files.")
    parser.add_argument('--sp_delimiter', '-d', default='_', help="Delimiter to split leaf names (default: '|').")
    parser.add_argument('--sp_field', '-f', type=int, default=0, help="Field index to extract after splitting (default: 0).")
    
    # Parse arguments
    args = parser.parse_args()

    # Iterate over input files
    for fname in args.files:
        # Read the tree from file
        t = Tree(open(fname))
        
        # Modify leaf names based on the delimiter and field index
        for leaf in t:
            leaf.name = leaf.name.split(args.sp_delimiter)[args.sp_field]
        
        # Print the modified tree
        print(t.write())

if __name__ == "__main__":
    main()
