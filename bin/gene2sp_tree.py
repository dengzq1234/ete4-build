#!/usr/bin/env python3

import sys
from ete4 import Tree

for fname in sys.argv[1:]:
    t = Tree(open(fname))
    for leaf in t:
        leaf.name = leaf.name.split('|')[0]

    print(t.write())