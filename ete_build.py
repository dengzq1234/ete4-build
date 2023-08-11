#!/usr/bin/env python

import nextflow

# pipeline = nextflow.Pipeline("main.nf")
# execution = pipeline.run(params={"param1": "123"})
# print(execution.status)

pipeline = nextflow.Pipeline('ete_build.nf')
execution = pipeline.run(params={
    "input": "../example.fasta",
    "output":"/home/deng/Projects/nextflow_projects/example_output/",
    "aligner": "mafft",
    "trimmer": "trimal",
    "tree_builder":"fasttree"
})

print(execution.status)
