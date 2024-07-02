#!/usr/bin/env nextflow

params.input = "$baseDir/data/"
params.output = "$baseDir/result"
params.thread = 1
params.aligner = "mafft"
params.trimmer = "trimal"
params.tree_builder = "fasttree"
params.memory = '4GB'
//memory_setting = params.memory ? params.memory : ''
params.time = '1h'
params.customConfig = null

bin = "$baseDir/bin"

// Load JSON config
def configPath = params.customConfig ?: "${baseDir}/config.json"
def config = file(configPath).text
def jsonConfig = new groovy.json.JsonSlurper().parseText(config)

if (file(params.input).isDirectory()) {
    FASTA_files = Channel.fromPath("${params.input}/*.{fa,faa,fasta}")
} else {
    FASTA_files = Channel.fromPath(params.input)
}

output_dir_structure = { fasta_name -> "${params.output}/${fasta_name}-${params.aligner}-${params.trimmer}-${params.tree_builder}" }

// Function to get MAFFT options
def getMafftOptions(alignConfig) {
    def methodConfig = alignConfig.methods[alignConfig.mode]
    def flag = methodConfig?.flag ?: "--auto"
    def options = "${flag} --ep ${alignConfig.ep} --op ${alignConfig.op} --maxiterate ${alignConfig.maxiterate}"
    
    // Add matrix options
    if (alignConfig.matrix == "BLOSUM") {
        options += " --bl ${alignConfig.blosum_coefficient}"
    } else if (alignConfig.matrix == "PAM") {
        options += " --jtt ${alignConfig.pam_coefficient}"
    }

    return options
}

// Function to get MUSCLE options
def getMuscleOptions(alignConfig) {
    def options = ""
    options += alignConfig.maxiters ? "-maxiters ${alignConfig.maxiters} " : ""
    options += alignConfig.diags ? "-diags " : ""
    return options
}

// Function to get T-Coffee options
def getTcoffeeOptions(alignConfig) {
    def options = "-n_core=${params.thread}"
    return options
}

// Function to get Clustal Omega options
def getClustaloOptions(alignConfig) {
    def options = ""
    options += alignConfig.dealign ? " --dealign" : ""
    options += alignConfig.mode == "full" ? " --full" : ""
    options += " --iterations ${alignConfig.iterations}"
    return options
}

// Function to get FAMSA options
def getFamsaOptions(alignConfig) {
    def options = ""
    options += " -gt ${alignConfig.gt}"
    options += " -r ${alignConfig.iterations}"
    return options
}

// Function to get Trimal options
def getTrimalOptions(trimConfig) {
    def options = ""
    options += trimConfig.gt ? "-gt ${trimConfig.gt} " : ""
    options += trimConfig.st ? "-st ${trimConfig.st} " : ""
    options += trimConfig.ct ? "-ct ${trimConfig.ct} " : ""
    return options
}

// Function to get Clipkit options
def getClipkitOptions(trimConfig) {
    def options = "--mode ${trimConfig.mode} "
    options += "--gaps ${trimConfig['gap-threshold']} "
    options += trimConfig.codon ? "--codon " : ""
    return options
}

// Function to get FastTree options
def getFastTreeOptions(buildConfig) {
    return ""
}

// Function to get PhyML options
// Function to get PhyML options
def getPhymlOptions(buildConfig) {
    def options = ""
    options += buildConfig.model ? "-m ${buildConfig.model} " : ""
    options += buildConfig.datatype ? "-d ${buildConfig.datatype} " : ""
    options += buildConfig.no_memory_check ? "--no_memory_check " : ""
    
    if (buildConfig.branch_support) {
        switch(buildConfig.branch_support) {
            case "bootstrap":
                options += "-b ${buildConfig.bootstrap_rep ?: 100} "
                if (buildConfig.bootstrap == "tbe") {
                    options += "--tbe "
                }
                break
            case "None":
                options += "-b 0 "
                break
            case "aLRT":
                options += "-b -1 "
                break
            case "Chi2":
                options += "-b -2 "
                break
            case "SH":
                options += "-b -4 "
                break
            case "Bayes":
                options += "-b -5 "
                break
            default:
                throw new Exception("Invalid branch support type: ${buildConfig.branch_support}")
        }
    }

    if (buildConfig.equilibrium_freq) {
        switch(buildConfig.equilibrium_freq) {
            case "empirical":
                options += "-f e "
                break
            case "ml":
                options += "-f m "
                break
            default:
                throw new Exception("Invalid equilibrium frequency: ${buildConfig.equilibrium_freq}")
        }
    }

    if (buildConfig.prop_invar) {
        if (buildConfig.prop_invar == "e") {
            options += "--pinv e "
        } else {
            options += "--pinv ${buildConfig.prop_invar} "
        }
    }

    if (buildConfig.gamma) {
        if (buildConfig.gamma == "e") {
            options += "--gamma e "
        } else {
            options += "--gamma ${buildConfig.gamma} "
        }
    }
    
    return options
}

// Function to get RAxML options
def getRaxmlOptions(buildConfig) {
    def options = ""
    options += buildConfig.algorithm ? "-f ${buildConfig.algorithm} " : ""
    options += buildConfig.r_seed ? "-p ${buildConfig.r_seed} " : ""
    options += buildConfig.model ? "-m ${buildConfig.model} " : ""
    options += buildConfig.bootstrap ? "-b ${buildConfig.bootstrap} " : ""
    return options
}

// Function to get IQ-TREE options
def getIqtreeOptions(buildConfig) {
    def options = ""
    options += buildConfig.st ? "-st ${buildConfig.st} " : ""
    options += buildConfig.alrt ? "-alrt ${buildConfig.alrt} " : ""
    options += buildConfig.seed ? "-seed ${buildConfig.seed} " : ""
    options += buildConfig.mode ? "-m ${buildConfig.mode} " : ""
    options += buildConfig.bootstrap_rep ? "-B ${buildConfig.bootstrap_rep} " : ""
    if (buildConfig.bootstrap == "tbe") {
        options += "--tbe "
    }
    return options
}


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
    echo $fasta_file
    awk '/^>/{if (seq) {print seq; seq=""} print \$0} {seq=seq\$0} END {print seq}' $fasta_file | \
    awk 'NR%2==0' | \
    awk '{if (length(\$0) > max) max = length(\$0); sum+=length(\$0)} END {print "Longest sequence: " max " characters"; print "Total sequences: " NR; print "Average sequence length: " sum/NR}'

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
    memory params.memory 
    time params.time
    errorStrategy 'retry'
    maxRetries 2
    publishDir path: { "${params.output}/${fasta_name}-${params.aligner}-${params.trimmer}-${params.tree_builder}" }, mode: 'copy'

    input:
    path fasta_file 

    output:
    path "${fasta_name}.aln.faa", emit: aln_seqs
    path "*.*", emit: aln_files
    stdout emit: align_stdout
    path "align.err", emit: align_err

    script:
    fasta_name = fasta_file.baseName
    def alignConfig = jsonConfig.aligner[params.aligner]
    
    def alignCmd = "mafft"
    def alignOptions = ""
    switch(params.aligner) {
        case "mafft":
            alignOptions = getMafftOptions(alignConfig)
            break
        case "muscle":
            alignCmd = "muscle"
            alignOptions = getMuscleOptions(alignConfig)
            break
        case "tcoffee":
            alignCmd = "t_coffee"
            alignOptions = getTcoffeeOptions(alignConfig)
            break
        case "clustalo":
            alignCmd = "clustalo"
            alignOptions = getClustaloOptions(alignConfig)
            break
        case "famsa":
            alignCmd = "famsa"
            alignOptions = getFamsaOptions(alignConfig)
            break
        default:
            throw new Exception("Invalid aligner: ${params.aligner}")
    }

    """
    num_sequences=\$(grep -c '^>' $fasta_file)
    start_time=\$(date +%s)
    echo "run ${alignCmd} with options: $alignOptions"
    if [ "${params.aligner}" == "mafft" ]; then
        ${alignCmd} ${alignOptions} --thread ${params.thread} $fasta_file > ${fasta_name}.aln.faa 2> align.err
    elif [ "${params.aligner}" == "muscle" ]; then
        ${alignCmd} ${alignOptions} -in $fasta_file -out ${fasta_name}.aln.faa 2> align.err
    elif [ "${params.aligner}" == "tcoffee" ]; then
        ${alignCmd} ${alignOptions} -in $fasta_file -outfile=${fasta_name}.aln.faa 2> align.err
    elif [ "${params.aligner}" == "clustalo" ]; then
        ${alignCmd} ${alignOptions} --thread ${params.thread} -i $fasta_file -o ${fasta_name}.aln.faa 2> align.err
    elif [ "${params.aligner}" == "famsa" ]; then
        ${alignCmd} ${alignOptions} -t ${params.thread} $fasta_file ${fasta_name}.aln.faa 2> align.err
    fi
    end_time=\$(date +%s)
    echo "Alignment $fasta_file took \$((end_time - start_time)) seconds."
    """
}

process trim {
    cpus params.thread
    memory params.memory
    time params.time
    errorStrategy 'retry'
    maxRetries 2
    publishDir path: { "${params.output}/${fasta_name}-${params.aligner}-${params.trimmer}-${params.tree_builder}" }, mode: 'copy'

    input:
    path aln_file 

    output:
    path "${fasta_name}.clean.alg.faa", emit: clean_aln_seqs
    path "*.*", emit: trim_files
    path "trim.out", emit: trim_out
    path "trim.err", emit: trim_err
    stdout emit: trim_stdout

    script:
    fasta_name = aln_file.baseName.replace(".aln", "")
    def trimConfig = jsonConfig.trimmer[params.trimmer]
    
    def trimCmd = "trimal"
    def trimOptions = ""
    switch(params.trimmer) {
        case "trimal":
            trimOptions = getTrimalOptions(trimConfig)
            break
        case "clipkit":
            trimCmd = "clipkit"
            trimOptions = getClipkitOptions(trimConfig)
            break
        default:
            throw new Exception("Invalid trimmer: ${params.trimmer}")
    }

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
    memory params.memory
    time params.time
    errorStrategy 'retry'
    maxRetries 2
    publishDir path: { "${params.output}/${fasta_name}-${params.aligner}-${params.trimmer}-${params.tree_builder}" }, mode: 'copy'
    
    input:
    path clean_aln_file 

    output:
    path "${fasta_name}.output.tree", emit: output_tree
    path "*.*", emit: build_files
    path "build.err", emit: build_err
    stdout emit: build_stdout

    script:
    fasta_name = clean_aln_file.baseName.replace(".clean.alg", "")
    def buildConfig = jsonConfig.tree_builder[params.tree_builder]
    def buildCmd = "FastTree"
    def buildOptions = ""
    switch(params.tree_builder) {
        case "fasttree":
            buildOptions = getFastTreeOptions(buildConfig)
            break
        case "phyml":
            buildCmd = "phyml"
            buildOptions = getPhymlOptions(buildConfig)
            break
        case "raxml":
            buildCmd = "raxmlHPC"
            buildOptions = getRaxmlOptions(buildConfig)
            break
        case "iqtree":
            buildCmd = "iqtree2"
            buildOptions = getIqtreeOptions(buildConfig)
            break
        default:
            throw new Exception("Invalid tree builder: ${params.tree_builder}")
    }

    """
    start_time=\$(date +%s)
    echo "run ${buildCmd} with options: $buildOptions"
    if [ "${params.tree_builder}" == "fasttree" ]; then
        ${buildCmd} ${buildOptions} $clean_aln_file > ${fasta_name}.output.tree 2> build.err
    elif [ "${params.tree_builder}" == "phyml" ]; then
        python ${bin}/FastaToPhylip.py $clean_aln_file && \
        ${buildCmd} ${buildOptions} -i ${fasta_name}.faa.clean.alg.phylip 2> build.err && \
        mv ${fasta_name}.faa.clean.alg.phylip_phyml_tree.txt ${fasta_name}.output.tree
    elif [ "${params.tree_builder}" == "raxml" ]; then
        ${buildCmd} ${buildOptions} -s $clean_aln_file -n ${fasta_name}.output.tree -T ${params.thread} 2> build.err
        cp RAxML_bestTree.${fasta_name}.output.tree ${fasta_name}.output.tree
    elif [ "${params.tree_builder}" == "iqtree" ]; then
        ${buildCmd} ${buildOptions} -s $clean_aln_file -T ${params.thread} 2> build.err
        cp ${fasta_name}.clean.alg.faa.treefile ${fasta_name}.output.tree
    fi
    end_time=\$(date +%s)
    echo "Tree building $fasta_name took \$((end_time - start_time)) seconds."
    """
}


workflow {
    FASTA_files.set{ fasta_ch }
    align(fasta_ch)
    trim(align.out.aln_seqs)
    build(trim.out.clean_aln_seqs)
    align.out.align_stdout.view { it -> println("[align] ${it}") }
    trim.out.trim_stdout.view { it -> println("[trim] ${it}") }
    build.out.build_stdout.view { it -> println("[build] ${it}") }
}

// workflow.onComplete {
//     file(".nextflow.log").moveTo("${params.output}/.nextflow.log")
// }
