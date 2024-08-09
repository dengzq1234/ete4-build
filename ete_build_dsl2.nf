#!/usr/bin/env nextflow
import groovy.json.JsonOutput

params.input = "$baseDir/data/"
params.output = "$baseDir/result"
params.thread = 1
params.aligner = "none" // "mafft"
params.trimmer = "none" // "trimal"
params.tree_builder = "none" // "fasttree"
params.memory = '4GB'
params.time = '1h'
params.customConfig = null

bin = "$baseDir/bin"

// Default configuration

def defaultConfig = [
    aligner: [
        mafft: [
            name: "mafft",
            mode: "auto",
            methods: [
                auto: [flag: "--auto"]
            ]
        ],
        muscle: [
            name: "muscle",
        ],
        tcoffee: [
            name: "t_coffee"
        ],
        clustalo: [
            name: "clustalo",
            dealign: false,
        ],
        famsa: [
            name: "famsa",
        ]
    ],
    trimmer: [
        trimal: [
            name: "trimal",
            gt: 0.1
        ],
        clipkit: [
            name: "clipkit",
            mode: "smart-gap",
            gap_threshold: 0.9,
            codon: false
        ],
        trim_alg_v2: [
            name: "trim_alg_v2.py",
            min_res_abs: 3,
            min_res_percent: 0.1
        ]
    ],
    tree_builder: [
        fasttree: [
            name: "fasttree"
        ],
        phyml: [
            name: "phyml",
            datatype: "aa",
            aa_models: ["LG", "WAG", "JTT", "MtREV", "Dayhoff", "DCMut", "RtREV", "CpREV", "VT", "AB", "Blosum62", "MtMam", "MtArt", "HIVw", "HIVb"],
            nt_models: ["HKY85", "JC69", "K80", "F81", "F84", "TN93", "GTR"],
            model: "JTT",
            no_memory_check: true,
            branch_support: "Chi2",
            equilibrium_freq: "empirical",
            prop_invar: "e",
            gamma: "e"
        ],
        raxml: [
            name: "raxmlHPC",
            algorithm: "d",
            r_seed: 31416,
            model: "PROTGAMMAJTT"
        ],
        iqtree: [
            name: "iqtree",
            st: "AA",
            alrt: 1000,
            seed: 31416,
            mode: "TESTONLY",
            bootstrap_rep: 1000,
            bootstrap: "fbp"
        ]
    ]
]

// Function to perform a deep copy of a map
def deepCopy(map) {
    new groovy.json.JsonSlurper().parseText(new groovy.json.JsonBuilder(map).toString())
}

// Function to load and merge custom config
def loadAndMergeConfig(defaultConfig, customConfigFile) {
    def config = deepCopy(defaultConfig)
    if (customConfigFile) {
        //println "Loading custom configuration from: ${customConfigFile}"
        def customConfig = new groovy.json.JsonSlurper().parseText(file(customConfigFile).text)
        //println "Custom Config: ${customConfig}"
        config = config + customConfig
    }
    println "workflow Config: ${config}"
    return config
}

// Ensure the configuration is loaded before any processes run
def jsonConfig = loadAndMergeConfig(defaultConfig, params.customConfig)

// Load custom config if provided old
// def jsonConfig = defaultConfig
// if (params.customConfig) {
//     def customConfig = new groovy.json.JsonSlurper().parseText(file(params.customConfig).text)
//     jsonConfig = defaultConfig + customConfig
// }

// Handle input files or directory
FASTA_files = file(params.input).isDirectory() ? Channel.fromPath("${params.input}/*.{fa,faa,fasta}") : Channel.fromPath(params.input)

output_dir_structure = { fasta_name -> "${params.output}/${fasta_name}-${params.aligner}-${params.trimmer}-${params.tree_builder}" }

// Function to get MAFFT options
def getMafftOptions(alignConfig) {
    def methodConfig = alignConfig.methods[alignConfig.mode]
    def flag = methodConfig?.flag ?: "--auto"
    def options = "${flag}"

    if (alignConfig.ep != null) {
        options += " --ep ${alignConfig.ep}"
    }
    if (alignConfig.op != null) {
        options += " --op ${alignConfig.op}"
    }
    if (alignConfig.maxiterate != null) {
        options += " --maxiterate ${alignConfig.maxiterate}"
    }

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
    options += alignConfig.iterations ? " --iterations ${alignConfig.iterations}": ""
    return options
}

// Function to get FAMSA options
def getFamsaOptions(alignConfig) {
    def options = ""
    options += alignConfig.gt ? " -gt ${alignConfig.gt}" : ""
    options += alignConfig.iterations ? " -r ${alignConfig.iterations}" : ""
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

def getTrimAlgV2Options(trimConfig) {
    def options = ""
    options += trimConfig.min_res_abs ? "--min_res_abs ${trimConfig.min_res_abs} " : ""
    options += trimConfig.min_res_percent ? "--min_res_percent ${trimConfig.min_res_percent} " : ""
    return options
}

// Function to get FastTree options
def getFastTreeOptions(buildConfig) {
    def options = ""

    // Handle datatype
    if (buildConfig.datatype == "nt") {
        options += " -nt"
    }
    // Handle model for aa
    if (buildConfig.datatype == "aa") {
        switch (buildConfig.model) {
            case "LG":
                options += " -lg"
                break
            case "WAG":
                options += " -wag"
                break
            // No need to add anything for JTT
        }
    }
    // Handle model for nt
    if (buildConfig.datatype == "nt") {
        switch (buildConfig.model) {
            case "GTR":
                options += " -gtr"
                break
            // No need to add anything for JC
        }
    }
    // Handle gamma
    if (buildConfig.gamma) {
        options += " -gamma"
    }
    // Handle pseudo
    if (buildConfig.pseudo != null) {
        options += " -pseudo ${buildConfig.pseudo}"
    }
    // Handle branch support
    if (buildConfig.branch_support) {
        options += " -boot ${buildConfig.bootstrap_rep ?: 1000}"
    }

    return options
}

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
    // println "IQ-TREE Options: ${buildConfig.bootstrap_rep}"
    if (buildConfig.bootstrap == "tbe") {
        options += "--tbe "
    }
    return options
}

// Function to filter only the used method in aligner
def filterUsedAlignerConfig(alignConfig) {
    def usedAlignerConfig = alignConfig.clone()
    if (alignConfig.methods) {
        usedAlignerConfig.methods = [(alignConfig.mode): alignConfig.methods[alignConfig.mode]]
    }
    return usedAlignerConfig
}

// Function to filter only the used options in trimmer
def filterUsedTrimmerConfig(trimConfig) {
    def usedTrimmerConfig = trimConfig.clone()
    return usedTrimmerConfig
}

// Function to filter only the used options in tree builder
def filterUsedTreeBuilderConfig(buildConfig) {
    def usedTreeBuilderConfig = buildConfig.clone()
    return usedTreeBuilderConfig
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
    if (params.aligner == "none") {
        // If no aligner is specified, just copy the input file to the output
        """
        cp $fasta_file ${fasta_name}.aln.faa
        """
    } else {
        def alignConfig = jsonConfig.aligner[params.aligner]
        
        
        if (!alignConfig) {
            throw new Exception("Aligner configuration for '${params.aligner}' is not found.")
        }

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

        // Additional debug print for alignOptions
        // println "Align Options: ${alignOptions}"
        """
        num_sequences=\$(grep -c '^>' $fasta_file)
        start_time=\$(date +%s)
        echo "run ${alignCmd} with options: $alignOptions"
        if [ "${params.aligner}" == "mafft" ]; then
            ${alignCmd} ${alignOptions} --thread ${params.thread} $fasta_file > ${fasta_name}.aln.faa 2> align.err
        elif [ "${params.aligner}" == "muscle" ]; then
            ${alignCmd} ${alignOptions} -align $fasta_file -output ${fasta_name}.aln.faa 2> align.err
        elif [ "${params.aligner}" == "tcoffee" ]; then
            ${alignCmd} ${alignOptions} -in $fasta_file -outfile=${fasta_name}.aln.faa 2> align.err
        elif [ "${params.aligner}" == "clustalo" ]; then
            ${alignCmd} ${alignOptions} --threads ${params.thread} -i $fasta_file -o ${fasta_name}.aln.faa 2> align.err
        elif [ "${params.aligner}" == "famsa" ]; then
            ${alignCmd} ${alignOptions} -t ${params.thread} $fasta_file ${fasta_name}.aln.faa 2> align.err
        fi
        end_time=\$(date +%s)
        echo "Alignment $fasta_file took \$((end_time - start_time)) seconds."
        """
    }
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
    path "*.*", emit: trim_files, optional: true
    stdout emit: trim_stdout
    path "trim.out",  optional: true
    path "trim.err", optional: true
    

    script:
    fasta_name = aln_file.baseName.replace(".aln", "")
    if (params.trimmer == "none") {
        // If no trimmer is specified, just copy the input file to the output
        println "No trimmer specified, copying the input file to the output."
        """
        cp $aln_file ${fasta_name}.clean.alg.faa
        """
    } else {
        def trimConfig = jsonConfig.trimmer[params.trimmer]
        def trimCmd = ""
        def trimOptions = ""
        switch(params.trimmer) {
            case "trimal":
                trimCmd = "trimal"
                trimOptions = getTrimalOptions(trimConfig)
                break
            case "clipkit":
                trimCmd = "clipkit"
                trimOptions = getClipkitOptions(trimConfig)
                break
            case "trim_alg_v2":
                trimCmd = "trim_alg_v2.py"
                trimOptions = getTrimAlgV2Options(trimConfig)
                break
            default:
                throw new Exception("Invalid trimmer: ${params.trimmer}")
        }

        """
        start_time=\$(date +%s)
        if [ "${params.trimmer}" == "trimal" ]; then
            echo "Running trimal!"
            ${trimCmd} ${trimOptions} -in $aln_file -out ${fasta_name}.clean.alg.faa -fasta 1> trim.out 2> trim.err
        elif [ "${params.trimmer}" == "clipkit" ]; then
            echo "Running clipkit!"
            ${trimCmd} ${trimOptions} $aln_file -o ${fasta_name}.clean.alg.faa 1> trim.out 2> trim.err
        elif [ "${params.trimmer}" == "trim_alg_v2" ]; then
            echo "Running trim_alg_v2!"
            python ${bin}/${trimCmd} ${trimOptions} -i $aln_file -o ${fasta_name}.clean.alg.faa 1> trim.out 2> trim.err
        fi
        end_time=\$(date +%s)
        echo "Trimming $fasta_name took \$((end_time - start_time)) seconds."
        """
    }
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
    if (params.tree_builder == "none") {
        // If no tree builder is specified, just copy the input file to the output
        """
        cp $clean_aln_file ${fasta_name}.output.tree
        """
    } else {
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
            ${buildCmd} ${buildOptions} -i ${fasta_name}.clean.alg.phylip 2> build.err && \
            mv ${fasta_name}.clean.alg.phylip_phyml_tree.txt ${fasta_name}.output.tree
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
}

process listInputFiles {
    input:
    path fasta_file

    output:
    path fasta_file

    script:
    """
    echo "Processing file: $fasta_file"
    """
}

process outputUsedConfig {
    input:
    val alignConfig
    val trimConfig
    val buildConfig

    output:
    path "used_config.json"

    publishDir params.output, mode: 'copy'

    script:
    // Apply filtering to the configurations
    def filteredAlignConfig = alignConfig ? filterUsedAlignerConfig(alignConfig) : null
    def filteredTrimConfig = trimConfig ? filterUsedTrimmerConfig(trimConfig) : null
    def filteredBuildConfig = buildConfig ? filterUsedTreeBuilderConfig(buildConfig) : null

    def usedConfig = [:]
    if (filteredAlignConfig) {
        usedConfig.aligner = [(params.aligner): filteredAlignConfig]
    }
    if (filteredTrimConfig) {
        usedConfig.trimmer = [(params.trimmer): filteredTrimConfig]
    }
    if (filteredBuildConfig) {
        usedConfig.tree_builder = [(params.tree_builder): filteredBuildConfig]
    }

    def jsonOutput = groovy.json.JsonOutput.toJson(usedConfig)
    """
    echo '${jsonOutput}' > used_config.json
    """
}


workflow {

    // FASTA_files | listInputFiles | view { it -> println("[listInputFiles] ${it}") }

    // Process all input files through the align, trim, and build processes
    FASTA_files
        .ifEmpty { error "No input files found in the provided input path: ${params.input}" }
        .set { parsed_files }
    
    // Retrieve the configurations used
    def alignConfig = params.aligner != "none" ? jsonConfig.aligner[params.aligner] : [:]
    def trimConfig = params.trimmer != "none" ? jsonConfig.trimmer[params.trimmer] : [:]
    def buildConfig = params.tree_builder != "none" ? jsonConfig.tree_builder[params.tree_builder] : [:]
    println "Align Config: ${alignConfig}"
    println "${params.trimmer}"
    println "Trim Config: ${trimConfig}"
    println "Build Config: ${buildConfig}"
    // Output the used configurations to a JSON file

    outputUsedConfig(alignConfig, trimConfig, buildConfig)

    align(parsed_files)
    trim(align.out.aln_seqs)
    build(trim.out.clean_aln_seqs)
    
    align.out.align_stdout.view { it -> println("[align] ${it}") }
    trim.out.trim_stdout.view { it -> println("[trim] ${it}") }
    build.out.build_stdout.view { it -> println("[build] ${it}") }
}

// workflow.onComplete {
//     file(".nextflow.log").moveTo("${params.output}/.nextflow.log")
// }
