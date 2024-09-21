import argparse

def create_mrbayes_commands(fasta_name, build_config, outfile, aln_type):
    """
    Generates MrBayes commands and writes them to a specified file.
    """
    # Generate MrBayes options formatted for the command string
    mb_options = f"""
    ngen={build_config['ngen']}
    nchains={build_config['nchains']}
    nruns={build_config['nruns']}
    samplefreq={build_config['samplefreq']}
    printfreq={build_config['printfreq']}
    diagnfreq={build_config['diagnfreq']}
    burninfrac={build_config['burninfrac']}
    stoprule={build_config['stoprule']}
    append={build_config['append']}
    """.strip().replace("\n", " ")  # Formatting options into a single line

    # Define the content of commands.txt using Python string formatting
    if aln_type == "DNA" or aln_type == "RNA":
        commands = f"""
    begin mrbayes;
        set autoclose=yes nowarn=yes;
        execute {fasta_name}.clean.alg.nex;
        lset nst=6 rates=invgamma;  # GTR model with gamma distribution and a proportion of invariable sites
        set seed={build_config['seed']} swapseed={build_config['swapseed']};
        mcmc {mb_options};
        sump burninfrac={build_config['burninfrac']};
        sumt burninfrac={build_config['burninfrac']};
    end;
        """
    else:
        commands = f"""
    begin mrbayes;
        set autoclose=yes nowarn=yes;
        execute {fasta_name}.clean.alg.nex;
        prset aamodelpr=fixed(wag);  # Choose the amino acid model here, e.g., wag
        set seed={build_config['seed']} swapseed={build_config['swapseed']};
        mcmc {mb_options};
        sump burninfrac={build_config['burninfrac']};
        sumt burninfrac={build_config['burninfrac']};
    end;
        """


    # Write the commands to the specified output file
    with open(outfile, 'w') as file:
        file.write(commands)

    print(f"MrBayes commands.txt created at {outfile} with the following content:")
    print(commands)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate MrBayes commands.txt")
    parser.add_argument("--fasta_name", required=True, help="Base name of the fasta file without extension")
    parser.add_argument("--outfile", required=True, help="Output file path for commands.txt")
    parser.add_argument("--aln_type", required=True, help="Alignment type: DNA, RNA, or protein")
    # MrBayes-specific arguments
    parser.add_argument("--ngen", type=int, required=True, help="Number of generations")
    parser.add_argument("--nchains", type=int, required=True, help="Number of chains")
    parser.add_argument("--nruns", type=int, required=True, help="Number of runs")
    parser.add_argument("--samplefreq", type=int, required=True, help="Sample frequency")
    parser.add_argument("--printfreq", type=int, required=True, help="Print frequency")
    parser.add_argument("--diagnfreq", type=int, required=True, help="Diagnostics frequency")
    parser.add_argument("--burninfrac", type=float, required=True, help="Burn-in fraction")
    parser.add_argument("--append", type=str, default='no', help="Append to the last checkpoint if true")
    parser.add_argument("--stoprule", type=str, default='no', help="Enable stop rule if true")
    parser.add_argument("--seed", type=int, required=True, help="Random seed for MrBayes")
    parser.add_argument("--swapseed", type=int, required=True, help="Swap seed for MrBayes")

    args = parser.parse_args()

    # Build the config dictionary from parsed arguments
    build_config = {
        'ngen': args.ngen,
        'nchains': args.nchains,
        'nruns': args.nruns,
        'samplefreq': args.samplefreq,
        'printfreq': args.printfreq,
        'diagnfreq': args.diagnfreq,
        'burninfrac': args.burninfrac,
        'append': args.append,
        'stoprule': args.stoprule,
        'seed': args.seed,
        'swapseed': args.swapseed,
    }

    # Create MrBayes commands file
    create_mrbayes_commands(args.fasta_name, build_config, args.outfile, args.aln_type)
