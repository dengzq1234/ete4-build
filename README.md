# Scalable Phylogenetic Tree Inference Nextflow Workflow Runner

This script provides a convenient way to run a Nextflow workflow with various configurations and parameters.

## Requirements

- Nextflow
- Python 3.x

## Usage

\```bash
python script_name.py --mode [local/slurm] --input INPUT_PATH --output OUTPUT_PATH [OPTIONS]
\```

### Arguments:

- `--mode`: Execution mode. Choose between `local` and `slurm`. This is a required argument.
- `--input`: Path to the input fasta file or directory. This is a required argument.
- `--output`: Path to the output directory where results will be saved. This is a required argument.
- `--aligner`: Alignment tool to use. Default is `mafft`.
- `--trimmer`: Trimming tool to use. Default is `trimal`.
- `--tree_builder`: Tree building tool to use. Default is `fasttree`.

### SLURM-specific arguments:

If you choose `slurm` as the execution mode, you need to provide additional SLURM-specific arguments:

- `--partition`: SLURM partition name. This is required if mode is `slurm`.
- `--time`: Time limit for SLURM jobs. Default is `1h`.
- `--memory`: Memory allocation for SLURM jobs. Default is `4GB`.
- `--cpus`: Number of CPUs for SLURM jobs. Default is `4`.

## Functionality

The script first generates a configuration file for Nextflow based on the provided execution mode and other parameters. It then runs the specified Nextflow workflow script (`ete_build_dsl2.nf` by default) with the given parameters.

The workflow processes the input fasta files to perform alignment, trimming, and tree building based on the specified tools.

## Output

The results will be saved in the specified output directory. The directory structure will be organized based on the input fasta file names and the tools used.