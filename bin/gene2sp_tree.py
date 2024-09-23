#!/usr/bin/env python3

import sys
from ete4 import Tree

sp_delimiter = '|'
sp_field = 0
for fname in sys.argv[1:]:
    t = Tree(open(fname))
    for leaf in t:
        leaf.name = leaf.name.split(sp_delimiter)[sp_field]

    print(t.write())