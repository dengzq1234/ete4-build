import argparse

def parse_nullable_int(value):
    """Convert 'null' or empty values to None; otherwise, parse as integer."""
    if value.lower() == 'null' or value == '':
        return None
    return int(value)

def parse_nullable_float(value):
    """Convert 'null' or empty values to None; otherwise, parse as float."""
    if value.lower() == 'null' or value == '':
        return None
    return float(value)

def parse_nullable_bool(value):
    """Convert 'null' or empty values to None; otherwise, parse as boolean."""
    if value.lower() == 'null' or value == '':
        return 'no'
    return 'yes' if value.lower() == 'true' else 'no'

def create_mrbayes_commands(fasta_name, build_config, outfile, aln_type):
    """
    Generates MrBayes commands and writes them to a specified file.
    """
    # Check if optional mcmc parameters are present, otherwise create an empty mcmc command
    mcmc_options = " ".join([
        f"ngen={build_config['ngen']}" if build_config['ngen'] is not None else "",
        f"nchains={build_config['nchains']}" if build_config['nchains'] is not None else "",
        f"nruns={build_config['nruns']}" if build_config['nruns'] is not None else "",
        f"samplefreq={build_config['samplefreq']}" if build_config['samplefreq'] is not None else "",
        f"printfreq={build_config['printfreq']}" if build_config['printfreq'] is not None else "",
        f"diagnfreq={build_config['diagnfreq']}" if build_config['diagnfreq'] is not None else "",
        f"burninfrac={build_config['burninfrac']}" if build_config['burninfrac'] is not None else "",
        f"stoprule={build_config['stoprule']}",
        f"append={build_config['append']}",
    ]).strip()

    # Generate lset options for DNA/RNA
    lset_options = "lset nst=6 rates=invgamma;" if aln_type in ["DNA", "RNA"] else ""

    # Generate prset options for proteins
    prset_options = "prset aamodelpr=fixed(wag);" if aln_type == "protein" else ""

    # Set command options
    set_options = f"set seed={build_config['seed']} swapseed={build_config['swapseed']};" if build_config['seed'] is not None and build_config['swapseed'] is not None else ""

    # Sump and Sumt options
    sump_options = f"sump burninfrac={build_config['burninfrac']};" if build_config['burninfrac'] is not None else "sump;"
    sumt_options = f"sumt burninfrac={build_config['burninfrac']};" if build_config['burninfrac'] is not None else "sumt;"

    # Define the full command based on options provided
    commands = f"""
    begin mrbayes;
        set autoclose=yes nowarn=yes;
        execute {fasta_name}.clean.alg.nex;
        {lset_options}
        {prset_options}
        {set_options}
        mcmc {mcmc_options};
        {sump_options}
        {sumt_options}
    end;
    """

    # Remove empty lines and extra spaces
    commands = "\n".join([line.strip() for line in commands.strip().splitlines() if line.strip()])

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
    parser.add_argument("--ngen", type=parse_nullable_int, default=None, help="Number of generations")
    parser.add_argument("--nchains", type=parse_nullable_int, default=None, help="Number of chains")
    parser.add_argument("--nruns", type=parse_nullable_int, default=None, help="Number of runs")
    parser.add_argument("--samplefreq", type=parse_nullable_int, default=None, help="Sample frequency")
    parser.add_argument("--printfreq", type=parse_nullable_int, default=None, help="Print frequency")
    parser.add_argument("--diagnfreq", type=parse_nullable_int, default=None, help="Diagnostics frequency")
    parser.add_argument("--burninfrac", type=parse_nullable_float, default=None, help="Burn-in fraction")
    parser.add_argument("--append", type=parse_nullable_bool, default='no', help="Append to the last checkpoint if true")
    parser.add_argument("--stoprule", type=parse_nullable_bool, default='no', help="Enable stop rule if true")
    parser.add_argument("--seed", type=parse_nullable_int, default=None, help="Random seed for MrBayes")
    parser.add_argument("--swapseed", type=parse_nullable_int, default=None, help="Swap seed for MrBayes")

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
