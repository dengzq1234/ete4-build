#!/usr/bin/env nextflow

params.input              = "$baseDir/data/NUP62.aa.fa"
params.output             = "$baseDir/result"
params.thread             = 4
params.aligner            = "mafft"
params.trimmer            = "trimal"
params.tree_builder       = "fasttree"
bin                       = "$baseDir/bin"

input_seqs   = file(params.input)
output_file  = params.output
thread       = params.thread 
aln_mode     = params.aligner
trim_mode    = params.trimmer 
builder_mode = params.tree_builder

process align {
    publishDir "${output_file}/${aln_mode}-${trim_mode}-${builder_mode}"

    input:
    path 'input.fa' 

    output:
    path 'input.aln.faa', emit: aln_seqs

    script:
    """
    if [[ "${aln_mode}" == 'mafft' ]]; then
        echo "run mafft!"
        mafft --anysymbol --thread ${thread} input.fa  > input.aln.faa
    elif [[ "${aln_mode}" == "tcoffee" ]]; then
        echo "run tcoffee!"
        t_coffee -cpu ${thread} -seq input.fa -output fasta_aln -outfile input.aln.faa
    elif [[ "${aln_mode}" == "clustalo" ]]; then
        echo "run clustalo!"
        clustalo --threads ${thread} --in input.fa -o input.aln.faa
    else
        echo "Invalid alignment mode: ${aln_mode}"
        exit 1
    fi
    """
}

process trim {
    publishDir "${output_file}/${aln_mode}-${trim_mode}-${builder_mode}"

    input:
    path x 

    output:
    path 'clean.alg.faa', emit: clean_aln_seqs

    script:
    """
    if [[ "${trim_mode}" == "trimal" ]]; then
        echo "trim!"
        trimal -in input.aln.faa -out clean.alg.faa -fasta -gt 0.1
    else
        cp input.aln.faa clean.alg.faa
    fi
    """
}

process build {
    publishDir "${output_file}/${aln_mode}-${trim_mode}-${builder_mode}"
    
    input:
    path x 

    output:
    path 'output.tree', emit: output_tree

    script:
    """
    if [[ "${builder_mode}" == "fasttree" ]]; then
        FastTree $x > output.tree
    elif [[ "${builder_mode}" == "phyml" ]]; then
        python ${bin}/FastaToPhylip.py $x
        phyml --model JTT --no_memory_check  --quiet  --pinv e --alpha e --nclasses 4 -o tlr -f m --bootstrap -2 --datatype aa --input clean.alg.phylip
        mv clean.alg.phylip_phyml_tree.txt output.tree
    elif [[ "${builder_mode}" == 'raxml' ]]; then
        raxmlHPC -f d -p 31416 -s $x -m PROTGAMMAJTT -n output.tree
        cp RAxML_bestTree.output.tree output.tree
    elif [[ "${builder_mode}" == 'iqtree' ]]; then
        iqtree -nt 1 -st AA -alrt 1000 -seed 31416 -s $x
        cp clean.alg.faa.treefile output.tree
    else
        echo "Invalid alignment mode: ${builder_mode}"
        exit 1
    fi
    """
}

// process organize_outputs {
//     output:
//     path 'final_output.txt'

//     script:
//     """
//     mkdir -p ${output_file}/${aln_mode}-${trim_mode}-${builder_mode}
//     mv result/${aln_mode}-${trim_mode}-${builder_mode}/* ${output_file}/${aln_mode}-${trim_mode}-${builder_mode}/
//     echo "Organized outputs"
//     touch final_output.txt
//     """
// }

workflow {
    align(input_seqs)
    trim(align.out.aln_seqs)
    build(trim.out.clean_aln_seqs)
    //organize_outputs()
}