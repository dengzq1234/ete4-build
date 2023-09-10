#!/usr/bin/env nextflow

params.input = "$baseDir/data/"  // Can be a directory or a single file
params.output = "$baseDir/result"
params.thread = 4
params.aligner = "mafft"
params.trimmer = "trimal"
params.tree_builder = "fasttree"
//params.memory = '4GB'
memory_setting = params.memory ? params.memory : ''
params.time = '1h'

bin = "$baseDir/bin"

// Check if input is a directory or a single file
if (file(params.input).isDirectory()) {
    FASTA_files = Channel.fromPath("${params.input}/*.fa")
} else {
    FASTA_files = Channel.fromPath(params.input)
}

// Define the output directory structure
output_dir_structure = { fasta_name -> "${params.output}/${fasta_name}-${params.aligner}-${params.trimmer}-${params.tree_builder}" }

process parseFasta {
    cpus 1
    memory '1GB'
    time '10m'
    errorStrategy 'retry'
    maxRetries 3
    publishDir path: { "${params.output}/${fasta_name}-${params.aligner}-${params.trimmer}-${params.tree_builder}" }, mode: 'copy'

    input:
    path fasta_file

    output:
        stdout emit: info

    script:
    fasta_name = fasta_file.baseName
    """
    # Using awk to parse the FASTA file
    echo $fasta_file
    awk '/^>/{if (seq) {print seq; seq=""} print \$0} {seq=seq\$0} END {print seq}' $fasta_file | \
    awk 'NR%2==0' | \
    awk '{if (length(\$0) > max) max = length(\$0); sum+=length(\$0)} END {print "Longest sequence: " max " characters"; print "Total sequences: " NR; print "Average sequence length: " sum/NR}'

    # Check for duplicated names
    grep '^>' $fasta_file | sort | uniq -d > duplicates.txt
    if [[ -s duplicates.txt ]]; then
        echo "Duplicated names found:"
        cat duplicates.txt
    else
        echo "No duplicated names found."
    fi
    rm duplicates.txt
    """
}

process align {
    cpus params.thread
    memory memory_setting
    time params.time
    errorStrategy 'retry'
    maxRetries 3
    publishDir path: { "${params.output}/${fasta_name}-${params.aligner}-${params.trimmer}-${params.tree_builder}" }, mode: 'copy'

    input:
    path fasta_file 

    output:
    path "${fasta_name}.aln.faa", emit: aln_seqs
    path "*.*", emit: aln_files  // Capture all files generated by the process
    //path "align.out", emit: align_out
    stdout emit: align_stdout
    path "align.err", emit: align_err

    script:
    fasta_name = fasta_file.baseName
    """
    start_time=\$(date +%s)
    if [[ "${params.aligner}" == 'mafft' ]]; then
        echo "run mafft!"
        mafft --anysymbol --thread ${params.thread} $fasta_file > ${fasta_name}.aln.faa 2> align.err
    elif [[ "${params.aligner}" == "tcoffee" ]]; then
        echo "run tcoffee!"
        t_coffee -cpu ${params.thread} -seq $fasta_file -output fasta_aln -outfile ${fasta_name}.aln.faa 2> align.err
    elif [[ "${params.aligner}" == "clustalo" ]]; then
        echo "run clustalo!"
        clustalo --threads ${params.thread} --in $fasta_file -o ${fasta_name}.aln.faa 2> align.err
    else
        echo "Invalid alignment mode: ${params.aligner}"
        exit 1
    fi
    end_time=\$(date +%s)
    echo "Alignment $fasta_file took \$((end_time - start_time)) seconds."
    """
}

process trim {
    cpus params.thread
    memory memory_setting
    time params.time
    errorStrategy 'retry'
    maxRetries 3
    publishDir path: { "${params.output}/${fasta_name}-${params.aligner}-${params.trimmer}-${params.tree_builder}" }, mode: 'copy'

    input:
    path aln_file 

    output:
    path "${fasta_name}.clean.alg.faa", emit: clean_aln_seqs
    path "*.*", emit: trim_files  // Capture all files generated by the process
    path "trim.out", emit: trim_out
    path "trim.err", emit: trim_err
    stdout emit: trim_stdout

    script:
    fasta_name = aln_file.baseName.replace(".aln", "")
    """
    start_time=\$(date +%s)
    if [[ "${params.trimmer}" == "trimal" ]]; then
        echo "trim!"
        trimal -in $aln_file -out ${fasta_name}.clean.alg.faa -fasta -gt 0.1 1> trim.out 2> trim.err
    else
        cp $aln_file ${fasta_name}.clean.alg.faa
    fi
    end_time=\$(date +%s)
    echo "Trimming $fasta_name took \$((end_time - start_time)) seconds."
    """
}

process build {
    cpus params.thread
    memory memory_setting
    time params.time
    errorStrategy 'retry'
    maxRetries 3
    publishDir path: { "${params.output}/${fasta_name}-${params.aligner}-${params.trimmer}-${params.tree_builder}" }, mode: 'copy'
    
    input:
    path clean_aln_file 

    output:
    path "${fasta_name}.output.tree", emit: output_tree
    path "*.*", emit: build_files  // Capture all files generated by the process
    //path "build.out", emit: build_out
    path "build.err", emit: build_err
    stdout emit: build_stdout

    script:
    fasta_name = clean_aln_file.baseName.replace(".clean.alg", "")
    """
    start_time=\$(date +%s)
    if [[ "${params.tree_builder}" == "fasttree" ]]; then
        FastTree $clean_aln_file > ${fasta_name}.output.tree  2> build.err
    elif [[ "${params.tree_builder}" == "phyml" ]]; then
        python ${bin}/FastaToPhylip.py $clean_aln_file
        phyml --model JTT --no_memory_check --quiet --pinv e --alpha e --nclasses 4 -o tlr -f m --bootstrap -2 --datatype aa --input ${fasta_name}.clean.alg.phylip -nt ${params.thread}  2> build.err
        mv ${fasta_name}.clean.alg.phylip_phyml_tree.txt ${fasta_name}.output.tree
    elif [[ "${params.tree_builder}" == 'raxml' ]]; then
        raxmlHPC -f d -p 31416 -s $clean_aln_file -m PROTGAMMAJTT -n ${fasta_name}.output.tree -T ${params.thread} 2> build.err
        cp RAxML_bestTree.${fasta_name}.output.tree ${fasta_name}.output.tree
    elif [[ "${params.tree_builder}" == 'iqtree' ]]; then
        iqtree -nt ${params.thread} -st AA -alrt 1000 -seed 31416 -s $clean_aln_file  2> build.err
        cp ${fasta_name}.clean.alg.faa.treefile ${fasta_name}.output.tree
    else
        echo "Invalid tree building mode: ${params.tree_builder}"
        exit 1
    fi
    end_time=\$(date +%s)
    echo "Tree building $fasta_name took \$((end_time - start_time)) seconds."
    """
}

workflow {
    FASTA_files.set{ fasta_ch }
    parseFasta(fasta_ch)
    align(fasta_ch)
    trim(align.out.aln_seqs)
    build(trim.out.clean_aln_seqs)
    parseFasta.out.info.view  { it -> println("[parseFasta] ${it}") }
    align.out.align_stdout.view { it -> println("[align] ${it}") }
    trim.out.trim_stdout.view { it -> println("[trim] ${it}") }
    build.out.build_stdout.view { it -> println("[build] ${it}") }
}

workflow.onComplete {
    file(".nextflow.log").moveTo("${params.output}/.nextflow.log")
}