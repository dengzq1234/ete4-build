#!/usr/bin/env nextflow

params.input              = "" //$baseDir/data/NUP62.aa.fa
params.output             = "" //$baseDir/result

params.thread             = 4
params.aligner            = "mafft" //tcoffee clustalo
params.trimmer            = "trimal"
params.tree_builder       = "iqtree" //raxml phyml

input_seqs   = file(params.input)
output_file  = params.output
thread       = params.thread 
aln_mode     = params.aligner
trim_mode    = params.trimmer 
builder_mode = params.tree_builder
bin          = "$baseDir/bin"

process align {
    publishDir "${output_file}/${aln_mode}-${trim_mode}-${builder_mode}"

    input:
    file 'input.fa' from input_seqs

    output:
    file 'input.aln.faa' into aln_seqs

    script:
    if (aln_mode == 'mafft')
        """
        echo "run mafft!"
        mafft --anysymbol --thread ${thread} input.fa  > input.aln.faa
        """
    else if (aln_mode == "tcoffee") 
        """
        echo "run tcoffee!"
        t_coffee -cpu ${thread} -seq input.fa -output fasta_aln -outfile input.aln.faa
        """
    else if (aln_mode == "clustalo")
        """
        echo "run clustalo!"
        clustalo --threads ${thread} --in input.fa -o input.aln.faa
        """
    else
        error "Invalid alignment mode: ${aln_mode}"
}

process trim {
    publishDir "${output_file}/${aln_mode}-${trim_mode}-${builder_mode}"

    input:
    file x from aln_seqs

    output:
    file 'clean.alg.faa' into clean_aln_seqs

    script:
    if (trim_mode == "trimal")
        """
        echo "trim!"
        trimal -in input.aln.faa -out clean.alg.faa -fasta -gt 0.1
        """
    else
        """
        cp input.aln.faa clean.alg.faa
        """

}

process build {
    publishDir "${output_file}/${aln_mode}-${trim_mode}-${builder_mode}"
    
    input:
    file x from clean_aln_seqs

    output:
    // stdout result
    file 'output.tree' into output_tree

    script:
    if (builder_mode == "fasttree")
        """
        FastTree $x > output.tree
        """
    else if (builder_mode == "phyml")
        """
        python ${bin}/FastaToPhylip.py $x
        phyml --model JTT --no_memory_check  --quiet  --pinv e --alpha e --nclasses 4 -o tlr -f m --bootstrap -2 --datatype aa --input clean.alg.phylip
        mv clean.alg.phylip_phyml_tree.txt output.tree
        """
    else if (builder_mode == 'raxml')
        """
        raxmlHPC -f d -p 31416 -s $x -m PROTGAMMAJTT -n output.tree
        cp RAxML_bestTree.output.tree output.tree
        """
    else if (builder_mode == 'iqtree')
        """
        iqtree -nt 1 -st AA -alrt 1000 -seed 31416 -s $x
        cp clean.alg.faa.treefile output.tree
        """
    else
        error "Invalid alignment mode: ${builder_mode}"

}

// process wrap {
//     input:
//     file input from clean_aln_seqs
//     file output from output_tree

//     output:
//     stdout result

//     script:
//     """
//     cp ${input} ${output_file}
//     cp ${output} ${output_file}
//     """

// }

// result.view { it.trim() }
//result.subscribe { println it }
