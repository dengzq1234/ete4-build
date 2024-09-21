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
params.supermatrix_mode = false // Whether to build a supermatrix
params.target_species = null // File path to target species list for supermatrix concatenation
bin = "$baseDir/bin"

// Default configuration

def defaultConfig = [
    aligner: [
        mafft: [
            name: "mafft",
            op: 1.53,
            ep: 0.123,
            maxiterate: 0
        ],
        muscle: [
            name: "muscle"
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
            gaps: 0.9,
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
            // datatype: 'aa',
            aa_model: "LG",  // Updated from the cfg
            nt_model: "HKY85",  // Updated from the cfg
            pinv: "e",  // Proportion of invariable sites, 'e' for estimation
            alpha: "e",  // Gamma distribution shape parameter, 'e' for estimation
            nclasses: 4,  // Number of rate categories
            optimisation: "tlr",  // Updated from the cfg
            frequencies: "m",  // Updated from the cfg
            bootstrap: -2,  // Updated to default Chi2-based parametric branch supports
            tbe: false,  // Disable TBE
            r_seed: 123456,  // Random seed
        ],
        raxml: [
            name: "raxmlHPC",
            algorithm: "d",
            r_seed: 31416,
            aa_model: "PROTGAMMAJTT",
            nt_model: "GTRGAMMA",
            bootstrap: 100,
        ],
        iqtree: [
            name: "iqtree",
            alrt: 1000,
            seed: 31416,
            model: "TESTONLY",
            tbe: false,  // Disable TBE
        ],
        mrbayes: [
            name: "mb",
            ngen: 100000,           // Number of generations
            nchains: 4,             // Number of chains
            nruns: 2,               // Number of runs
            nst: 1,                // Substitution model for dna
            rates: "equal",      // Rates model for dna Equal/Gamma/LNorm/Propinv/Invgamma/Adgamma/Kmixture 
            aamodelpr: "fixed(wag)", // Amino acid model
            diagnfreq: 5000,        // Frequency of diagnosing
            samplefreq: 500,        // Frequency of sampling
            printfreq: 1000,         // Frequency of printing
            burninfrac: 0.25,       // Burn-in fraction
            append: "no",           // Append to last checkpoint
            stoprule: "no",
            seed: 1726956368,                // Seed
            swapseed: 1726956368             // Swap seed
        ]
    ]
]

// Function to perform a deep copy of a map
def deepCopy(map) {
    new groovy.json.JsonSlurper().parseText(new groovy.json.JsonBuilder(map).toString())
}

// Function to load and merge custom config
// def loadAndMergeConfig(defaultConfig, customConfigFile) {
//     def config = [:]  // Start with an empty map

//     if (customConfigFile) {
//         // Load the custom configuration first
//         def customConfig = new groovy.json.JsonSlurper().parseText(file(customConfigFile).text)
//         config = config + customConfig
//         // println "check ${customConfig}"
//     } 
//     // Then add the default configuration, so customConfig has priority
//     config = defaultConfig + config
    
//     //println "workflow Config: ${config}"
//     return config
// }
// Function to load and merge custom config
def loadAndMergeConfig(defaultConfig, customConfigFile) {
    def config = [:]  // Start with an empty map

    if (customConfigFile) {
        // Load the custom configuration first
        def customConfig = new groovy.json.JsonSlurper().parseText(file(customConfigFile).text)

        // Only overwrite the aligner, trimmer, and tree_builder if they are not empty
        if (customConfig.containsKey('aligner') && !customConfig.aligner.isEmpty()) {
            config.aligner = customConfig.aligner
        } else {
            config.aligner = defaultConfig.aligner
        }

        if (customConfig.containsKey('trimmer') && !customConfig.trimmer.isEmpty()) {
            config.trimmer = customConfig.trimmer
        } else {
            config.trimmer = defaultConfig.trimmer
        }

        if (customConfig.containsKey('tree_builder') && !customConfig.tree_builder.isEmpty()) {
            config.tree_builder = customConfig.tree_builder
        } else {
            config.tree_builder = defaultConfig.tree_builder
        }

        // Handle other potential sections of the config that are not aligner, trimmer, or tree_builder
        customConfig.each { key, value ->
            if (!['aligner', 'trimmer', 'tree_builder'].contains(key)) {
                config[key] = value
            }
        }
    } else {
        // If no custom config is provided, just use the default
        config = defaultConfig
    }

    // Then add any other sections from the default config that aren't in the custom config
    config = defaultConfig + config

    // Return the final merged configuration
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
    def options = ""

    // Handle basic parameters
    if (alignConfig.ep != null) {
        options += " --ep ${alignConfig.ep}"
    }
    if (alignConfig.op != null) {
        options += " --op ${alignConfig.op}"
    }
    if (alignConfig.maxiterate != null) {
        options += " --maxiterate ${alignConfig.maxiterate}"
    }
    if (alignConfig.retree != null) {
        options += " --retree ${alignConfig.retree}"
    }

    // Add matrix options only if matrix is specified as "BLOSUM" or "PAM"
    if (alignConfig.matrix == "BLOSUM") {
        options += " --bl ${alignConfig.blosum_coefficient}"
    } else if (alignConfig.matrix == "PAM") {
        options += " --jtt ${alignConfig.pam_coefficient}"
    }

    // Handle boolean flags
    println "AlignConfig: ${alignConfig.auto}"
    if (alignConfig.auto) {
        options += " --auto"
    }
    if (alignConfig.localpair) {
        options += " --localpair"
    }
    if (alignConfig.globalpair) {
        options += " --globalpair"
    }
    if (alignConfig.genafpair) {
        options += " --genafpair"
    }
    if (alignConfig.nofft) {
        options += " --nofft"
    }
    if (alignConfig.parttree) {
        options += " --parttree"
    }

    println "MAFFT Options: ${options}"
    return options
}

// Function to get MUSCLE options
def getMuscleOptions(alignConfig) {
    def options = ""
    if (alignConfig.replicates) {
        options += " -replicates ${alignConfig.replicates}"
    }
    if (alignConfig.perturb) {
        options += " -perturb ${alignConfig.perturb}"
    }
    if (alignConfig.perm) {
        options += " -perm ${alignConfig.perm}"
    }
    if (alignConfig.consiters){
        options += " -consiters ${alignConfig.consiters}"
    }
    if (alignConfig.refineiters) {
        options += " -refineiters ${alignConfig.refineiters}"
    }
    
    if (alignConfig.stratified) {
        options += " -stratified"
    }
    if (alignConfig.diversified) {
        options += " -diversified"
    }

    println "MUSCLE Options: ${options}"
    return options
}

// Function to get T-Coffee options
def getTcoffeeOptions(alignConfig) {
    def options = "-n_core=${params.thread}"
    return options
}

// Function to get Clustal Omega options
// def getClustaloOptions(alignConfig) {
//     def options = ""
//     options += alignConfig.dealign ? " --dealign" : ""
//     options += alignConfig.mode == "full" ? " --full" : ""
//     options += alignConfig.iterations ? " --iterations ${alignConfig.iterations}": ""
//     return options
// }

def getClustaloOptions(alignConfig) {
    def options = ""
    if (alignConfig.dealign) {
        options += " --dealign"
    }
    if (alignConfig.full) {
        options += " --full"
    }
    if (alignConfig.full_iter) {
        options += " --full-iter"
    }
    if (alignConfig.iterations) {
        options += " --iterations ${alignConfig.iterations}"
    }
    if (alignConfig.max_guidetree_iterations) {
        options += " --max-guidetree-iterations ${alignConfig.max_guidetree_iterations}"
    }
    if (alignConfig.max_hmm_iterations) {
        options += " --max-hmm-iterations ${alignConfig.max_hmm_iterations}"
    }
    println "Clustal Omega Options: ${options}"
    return options
}

// Function to get FAMSA options
def getFamsaOptions(alignConfig) {
    def options = ""

    // Handle guide tree (gt)
    if (alignConfig.gt) {
        options += " -gt ${alignConfig.gt}"
    }

    // Handle medoid tree option
    if (alignConfig.medoidtree ) {
        options += " -medoidtree"
    }

    // Handle refine mode
    if (alignConfig.refine_mode) {
        options += " -refine_mode ${alignConfig.refine_mode}"
    }

    // Handle refinement iterations (r)
    if (alignConfig.r != null) {
        options += " -r ${alignConfig.r}"
    }

    // Handle gap penalties
    if (alignConfig.go != null) {
        options += " -go ${alignConfig.go}"
    }
    if (alignConfig.ge != null) {
        options += " -ge ${alignConfig.ge}"
    }
    if (alignConfig.tgo != null) {
        options += " -tgo ${alignConfig.tgo}"
    }
    if (alignConfig.tge != null) {
        options += " -tge ${alignConfig.tge}"
    }

    // Handle gap cost scaler terms
    if (alignConfig.gsd != null) {
        options += " -gsd ${alignConfig.gsd}"
    }
    if (alignConfig.gsl != null) {
        options += " -gsl ${alignConfig.gsl}"
    }

    // Handle disabling options (dgr, dgo, dsp)
    if (alignConfig.dgr != null && alignConfig.dgr) {
        options += " -dgr"
    }
    if (alignConfig.dgo != null && alignConfig.dgo) {
        options += " -dgo"
    }
    if (alignConfig.dsp != null && alignConfig.dsp) {
        options += " -dsp"
    }
    println "Famsa Options: ${options}"
    return options
}


// Function to get Trimal options
def getTrimalOptions(trimConfig) {
    def options = ""

    // Handle gap threshold (gt)
    if (trimConfig.gt != null) {
        options += "-gt ${trimConfig.gt} "
    }

    // Handle minimum average similarity threshold (st)
    if (trimConfig.st != null) {
        options += "-st ${trimConfig.st} "
    }

    // Handle minimum percentage of positions to conserve (ct)
    if (trimConfig.ct != null) {
        options += "-ct ${trimConfig.ct} "
    }

    // Handle gappyout option
    if (trimConfig.gappyout == true) {
        options += "-gappyout "
    }

    // Handle sliding window size (w)
    if (trimConfig.w != null) {
        options += "-w ${trimConfig.w} "
    }

    // Handle strictplus option (for NJ tree reconstruction)
    if (trimConfig.strictplus == true) {
        options += "-strictplus "
    }

    // Handle automated1 option (for ML tree reconstruction)
    if (trimConfig.automated1 == true) {
        options += "-automated1 "
    }
    println "Trimal Options: ${options}"
    // Return the assembled options string
    return options
}

// Function to get Clipkit options
def getClipkitOptions(trimConfig) {
    def options = ""

    // Ensure mode is specified, otherwise use a default
    if (trimConfig.mode) {
        options += "--mode ${trimConfig.mode} "
    } else {
        options += "--mode smart-gap " // Default mode
    }

    // Handle gaps threshold
    if (trimConfig.gaps != null) {
        options += "--gaps ${trimConfig.gaps} "
    } else {
        options += "--gaps 0.9 " // Default gap threshold
    }

    // Handle codon mode
    if (trimConfig.codon) {
        options += "--codon "
    }
    println "Clipkit Options: ${options}"
    return options
}

def getTrimAlgV2Options(trimConfig) {
    def options = ""
    options += trimConfig.min_res_abs ? "--min_res_abs ${trimConfig.min_res_abs} " : ""
    options += trimConfig.min_res_percent ? "--min_res_percent ${trimConfig.min_res_percent} " : ""
    println "TrimAlgV2 Options: ${options}"
    return options
}

// Function to get FastTree options
def getFastTreeOptions(buildConfig) {
    def options = ""

    // Handle model for amino acids
    if (buildConfig.aa_model) {
        switch (buildConfig.aa_model) {
            case "LG":
                options += " -lg"
                break
            case "WAG":
                options += " -wag"
                break
            case "JTT":
                break
            // No need to add anything for the default model (JTT+CAT)
        }
    }

    // Handle model for nucleotides
    if (buildConfig.nt_model) {
        switch (buildConfig.nt_model) {
            case "GTR":
                options += " -gtr"
                break
            case "JC":
                break
        }
    }

    // Handle gamma distribution
    if (buildConfig.gamma) {
        options += " -gamma"
    }

    // Handle pseudo-likelihood support values
    if (buildConfig.pseudo) {
        options += " -pseudo"
    }

    // Handle bootstrap or nosupport
    if (buildConfig.bootstrap == 0) {
        options += " -nosupport"
    } else if (buildConfig.bootstrap) {
        options += " -boot ${buildConfig.bootstrap ?: 1000}"
    }
    // Handle SPR rounds (minimum-evolution SPR moves)
    if (buildConfig.spr != null) {
        options += " -spr ${buildConfig.spr}"
    }

    // Handle ML model accuracy categories (MLACC)
    if (buildConfig.mlacc != null) {
        options += " -mlacc ${buildConfig.mlacc}"
    }

    // Handle slow NNI moves
    if (buildConfig.slownni) {
        options += " -slownni"
    }
    println "FastTree Options: ${options}"
    return options
}

// Function to get PhyML options
def getPhymlOptions(buildConfig, aln_type) {
    def options = ""
    // check point for datatype
    // Handle amino acid models
    if (aln_type == 'aa') {
        if (buildConfig.aa_model) {
            switch (buildConfig.aa_model) {
                case "LG":
                    options += " -m LG"
                    break
                case "WAG":
                    options += " -m WAG"
                    break
                case "JTT":
                    options += " -m JTT"
                    break
                case "MtREV":
                    options += " -m MtREV"
                    break
                case "Dayhoff":
                    options += " -m Dayhoff"
                    break
                case "DCMut":
                    options += " -m DCMut"
                    break
                case "RtREV":
                    options += " -m RtREV"
                    break
                case "CpREV":
                    options += " -m CpREV"
                    break
                case "VT":
                    options += " -m VT"
                    break
                case "AB":
                    options += " -m AB"
                    break
                case "Blosum62":
                    options += " -m Blosum62"
                    break
                case "MtMam":
                    options += " -m MtMam"
                    break
                case "MtArt":
                    options += " -m MtArt"
                    break
                case "HIVw":
                    options += " -m HIVw"
                    break
                case "HIVb":
                    options += " -m HIVb"
                    break
                case "custom":
                    options += " -m custom"
                    break
            }
        }
    } else {
        // Handle nucleotide models
        if (buildConfig.nt_model) {
            switch (buildConfig.nt_model) {
                case "HKY85":
                    options += " -m HKY85"
                    break
                case "JC69":
                    options += " -m JC69"
                    break
                case "K80":
                    options += " -m K80"
                    break
                case "F81":
                    options += " -m F81"
                    break
                case "F84":
                    options += " -m F84"
                    break
                case "TN93":
                    options += " -m TN93"
                    break
                case "GTR":
                    options += " -m GTR"
                    break
                case "custom":
                    options += " -m custom"
                    break
            }
        }
    }
    

    

    // Handle proportion of invariable sites (pinv)
    if (buildConfig.pinv != null) {
        options += buildConfig.pinv == "e" ? " --pinv e" : " --pinv ${buildConfig.pinv}"
    }

    // Handle gamma distribution shape parameter (alpha)
    if (buildConfig.alpha != null) {
        options += buildConfig.alpha == "e" ? " --alpha e" : " --alpha ${buildConfig.alpha}"
    }

    // Handle number of rate categories (nclasses)
    if (buildConfig.nclasses) {
        options += " --nclasses ${buildConfig.nclasses}"
    }

    // Handle tree optimisation (optimisation)
    if (buildConfig.optimisation) {
        options += " -o ${buildConfig.optimisation}"
    }

    // Handle frequencies (frequencies)
    if (buildConfig.frequencies) {
        switch (buildConfig.frequencies) {
            case "e":
                options += " -f e"
                break
            case "m":
                options += " -f m"
                break
            case "o":
                options += " -f o"
                break
            default:
                options += " -f ${buildConfig.frequencies}"
                break
        }
    }

    // Handle bootstrap or branch support (-b)
    if (buildConfig.bootstrap) {
        options += " -b ${buildConfig.bootstrap}"
    }

    // Handle TBE instead of FBP if tbe is True
    if (buildConfig.tbe) {
        options += " --tbe"
    }

    // Handle random seed (-r)
    if (buildConfig.r_seed) {
        options += " --r_seed ${buildConfig.r_seed}"
    }

    options += " --quiet"
    options += " --no_memory_check"

    return options
}


// Function to get RAxML options
// def getRaxmlOptions(buildConfig) {
//     def options = ""
//     options += buildConfig.algorithm ? "-f ${buildConfig.algorithm} " : ""
//     options += buildConfig.r_seed ? "-p ${buildConfig.r_seed} " : ""
//     options += buildConfig.model ? "-m ${buildConfig.model} " : ""
//     options += buildConfig.bootstrap ? "-b ${buildConfig.bootstrap} " : ""
//     return options
// }

def getRaxmlOptions(buildConfig, aln_type){
    def options = ""
    if (buildConfig.algorithm) {
        options += "-f ${buildConfig.algorithm} "
    }
    
    if (aln_type == 'aa'){
        if (buildConfig.aa_model) {
            options += "-m ${buildConfig.aa_model} "
        }
    } else {
        if (buildConfig.nt_model) {
            options += "-m ${buildConfig.nt_model} "
        }
    }
    println "RAxML Options: ${aln_type}, ${buildConfig.nt_model} "
    if (buildConfig.r_seed) {
        options += "-p ${buildConfig.r_seed} "
    }

    if (buildConfig.bootstrap) {
        options += "-N ${buildConfig.bootstrap} "
    }
    return options
}

// Function to get IQ-TREE options
def getIqtreeOptions(buildConfig) {
    def options = ""
    
    if (buildConfig.alrt) {
        options += "-alrt ${buildConfig.alrt} "
    }

    if (buildConfig.seed) {
        options += "-seed ${buildConfig.seed} "
    }

    if (buildConfig.mode) {
        options += "-m ${buildConfig.mode} "
    }
    
    if (buildConfig.st) {
        options += "-st ${buildConfig.st} "
    }

    if (buildConfig.ufboot) {
        options += "-B ${buildConfig.ufboot} "
    }

    // println "IQ-TREE Options: ${buildConfig.bootstrap_rep}"
    if (buildConfig.tbe) {
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

def detectAlignmentType(alignmentFile) {
    // Read the alignment file as text
    def lines = file(alignmentFile).readLines()

    // Concatenate all sequence lines (ignoring lines that look like headers or gaps)
    def sequenceData = lines.findAll { line -> !line.startsWith(">") && line.trim() != "" }
                             .join("")
                             .replace("-", "")  // Ignore gaps ("-")

    // Regular expressions for nucleotides and amino acids (case-insensitive with (?i))
    def nt_regex = /(?i)^[ACGTURYKMSWBDHVN]+$/    // IUPAC codes for nucleotides
    def aa_regex = /(?i)^[ACDEFGHIKLMNPQRSTVWYBXZ]+$/ // IUPAC codes for amino acids
    
    // Check if the sequence data matches nucleotide or amino acid patterns
    if (sequenceData ==~ nt_regex) {
        return "nt"
    } else if (sequenceData ==~ aa_regex) {
        return "aa"
    } else {
        throw new Exception("Cannot determine alignment type. Neither nucleotide nor amino acid patterns match.")
    }
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
            ${alignCmd} -align $fasta_file -output ${fasta_name}.aln.faa ${alignOptions} 2> align.err
        elif [ "${params.aligner}" == "tcoffee" ]; then
            ${alignCmd} ${alignOptions} -in $fasta_file -output=fasta_aln -outfile=${fasta_name}.aln.faa 2> align.err
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

process concatSupermatrix {
    cpus params.thread
    memory params.memory
    time params.time
    errorStrategy 'retry'
    maxRetries 2
    publishDir path: { "${params.output}/supermatrix-${params.aligner}-${params.trimmer}-${params.tree_builder}" }, mode: 'copy'

    input:
    //path clean_aln_files 
    path clean_aln_files // Collect all files into one input
    
    output:
    path "supermatrix.clean.alg.faa", emit: supermatrix_output
    path "*.*", emit: concat_files

    script:
    if (!params.target_species) {
        throw new Exception("Target species file is required when supermatrix mode is enabled.")
    }
    def targetSpeciesFile = file(params.target_species).toAbsolutePath().toString()

    // Debug: Print the collected alignment files for verification
    println "Collected alignment files for supermatrix concatenation:"
    clean_aln_files.each { println it.toString() }

    """
    python ${bin}/concat_aln.py -a ${clean_aln_files.join(" ")} --taxa ${targetSpeciesFile} -o supermatrix.clean.alg.faa -p partition_file.txt
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
    
    def aln_name = clean_aln_file.baseName.replace(".clean.alg", "")
    
    // Construct the correct file path for alignment
    def aln_file_path = file("${params.output}/${aln_name}-${params.aligner}-${params.trimmer}-${params.tree_builder}/${clean_aln_file.name}")
    
    
    def aln_type = detectAlignmentType(aln_file_path)
    def workDir = task.workDir.toString() // Get the current working directory for the task
    
    //println "Building tree for: ${clean_aln_file}, detected type: ${aln_type}"
    
    fasta_name = clean_aln_file.baseName.replace(".clean.alg", "") // for the bash script
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
                buildOptions = getPhymlOptions(buildConfig, aln_type)
                break
            case "raxml":
                buildCmd = "raxmlHPC"
                buildOptions = getRaxmlOptions(buildConfig, aln_type)
                break
            case "iqtree":
                buildCmd = "iqtree2"
                buildOptions = getIqtreeOptions(buildConfig)
                break
            case "mrbayes":
                buildCmd = "mb"
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
            # Convert the input FASTA file to PHYLIP format using the updated script
            python ${bin}/FastaToPhylip.py fasta2phylip -i $clean_aln_file -o ${fasta_name}.clean.alg.phylip && \
            # Run PhyML on the converted PHYLIP file
            ${buildCmd} ${buildOptions} -i ${fasta_name}.clean.alg.phylip 2> build.err && \
            # Move the output tree file to the desired output location
            mv ${fasta_name}.clean.alg.phylip_phyml_tree.txt ${fasta_name}.output.tree
        elif [ "${params.tree_builder}" == "raxml" ]; then
            ${buildCmd} ${buildOptions} -s $clean_aln_file -n ${fasta_name}.output.tree -T ${params.thread} 2> build.err
            cp RAxML_bestTree.${fasta_name}.output.tree ${fasta_name}.output.tree
        elif [ "${params.tree_builder}" == "iqtree" ]; then
            ${buildCmd} ${buildOptions} -s $clean_aln_file -T ${params.thread} 2> build.err
            cp ${fasta_name}.clean.alg.faa.treefile ${fasta_name}.output.tree
        elif [ "${params.tree_builder}" == "mrbayes" ]; then
            # Create the commands.txt file for MrBayes
            python ${bin}/FastaToPhylip.py fasta2nexus -i $clean_aln_file -o ${fasta_name}.clean.alg.nex 
            # Generate MrBayes commands.txt using the Python script
            python ${bin}/mrbayes_command.py \
                --fasta_name ${fasta_name} \
                --outfile commands.txt \
                --aln_type ${aln_type} \
                --ngen ${buildConfig.ngen} \
                --nchains ${buildConfig.nchains} \
                --nruns ${buildConfig.nruns} \
                --samplefreq ${buildConfig.samplefreq} \
                --printfreq ${buildConfig.printfreq} \
                --burninfrac ${buildConfig.burninfrac} \
                --diagnfreq ${buildConfig.diagnfreq} \
                --append ${buildConfig.append} \
                --seed ${buildConfig.seed} \
                --swapseed ${buildConfig.swapseed}

            # Run MrBayes using the generated commands.txt
            ${buildCmd} < commands.txt
            cp ${fasta_name}.clean.alg.nex.tre ${fasta_name}.output.tree  
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

process outputUsedConfigAsCfg {
    input:
    val alignConfig
    val trimConfig
    val buildConfig

    output:
    path "used_config.cfg"

    publishDir params.output, mode: 'copy'

    script:
    def configToCfg = { toolType, toolConfigs ->
        def cfgContent = ""
        toolConfigs.each { toolName, toolConfig ->
            def sectionName = "[${toolName}_default]"
            cfgContent += "${sectionName}\n"
            cfgContent += "_desc = '${toolName.capitalize()} with default parameters'\n"
            cfgContent += "_app = ${toolName}\n"
            
            toolConfig.each { key, value ->
                if (key == "name" || value instanceof Map) {
                    return
                }
                if (value instanceof Boolean) {
                    value = value ? "True" : "False"
                }
                cfgContent += "${key} = ${value}\n"
            }

            if (toolConfig.methods) {
                toolConfig.methods.each { method, methodSettings ->
                    sectionName = "[${toolName}_${method}]"
                    cfgContent += "\n${sectionName}\n"
                    cfgContent += "_inherits = ${toolName}_default\n"
                    cfgContent += "_desc = '${toolName.capitalize()} with ${method} method'\n"
                    methodSettings.each { methodKey, methodValue ->
                        if (methodValue instanceof Boolean) {
                            methodValue = methodValue ? "True" : "False"
                        }
                        cfgContent += "${methodKey} = ${methodValue}\n"
                    }
                }
            }

            cfgContent += "\n"
        }
        return cfgContent
    }

    def usedConfig = [:]
    if (alignConfig) {
        usedConfig.aligner = [(params.aligner): alignConfig]
    }
    if (trimConfig) {
        usedConfig.trimmer = [(params.trimmer): trimConfig]
    }
    if (buildConfig) {
        usedConfig.tree_builder = [(params.tree_builder): buildConfig]
    }

    def cfgOutput = ""
    usedConfig.each { toolType, toolConfigs ->
        cfgOutput += configToCfg(toolType, toolConfigs)
    }

    """
    echo '${cfgOutput}' > used_config.cfg
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

    outputUsedConfig(alignConfig, trimConfig, buildConfig)
    outputUsedConfigAsCfg(alignConfig, trimConfig, buildConfig)

    // process to to either gene tree or speceis tree
    align(parsed_files)
    trim(align.out.aln_seqs)
    
    if (params.supermatrix_mode) {
        // Collect all trimmed alignment files into a single set
        def collected_trimmed_files = trim.out.clean_aln_seqs.collect()

        // Pass the collected files to concatSupermatrix
        concatSupermatrix(collected_trimmed_files)

        // Use the resulting supermatrix output for tree building
        build(concatSupermatrix.out.supermatrix_output)
    } else {
        // Build trees individually from each trimmed alignment
        build(trim.out.clean_aln_seqs)
    }
    
    align.out.align_stdout.view { it -> println("[align] ${it}") }
    trim.out.trim_stdout.view { it -> println("[trim] ${it}") }
    build.out.build_stdout.view { it -> println("[build] ${it}") }
}

// workflow.onComplete {
//     file(".nextflow.log").moveTo("${params.output}/.nextflow.log")
// }
