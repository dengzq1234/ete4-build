#!/usr/bin/env python

import argparse
import os
import sys
#import nextflow
import subprocess

# Predefined workflows
PREDEFINED_WORKFLOWS = {
    "workflow1": {"aligner": "mafft", "trimmer": "trimal", "tree_builder": "fasttree"},
    "workflow2": {"aligner": "famsa", "trimmer": "trimal", "tree_builder": "fasttree"},
    "ana-workflow": {"aligner": "hybrid", "trimmer": "trimal", "tree_builder": "fasttree"},
    # Add more predefined workflows as needed
}

def generate_nextflow_config(execution_mode, partition=None, time=None, memory=None, cpus=None):
    if execution_mode == "slurm":
        config_content = f"""
        profiles {{
            slurm {{
                process {{
                    executor = 'slurm'
                    queue = '{partition}'
                    time = '{time}'
                    memory = '{memory}'
                    cpus = {cpus}
                }}
            }}
        }}
"""
    elif execution_mode == "local":
        config_content = """
profiles {
    local {
        process {
            executor = 'local'
        }
    }
}
"""
    else:
        raise ValueError(f"Unsupported execution mode: {execution_mode}")

    with open("nextflow.config", "w") as f:
        f.write(config_content)

def run_nextflow(mode, input_file, output_dir, aligner, trimmer, tree_builder, resume=False, script="ete_build_dsl2.nf"):
    generate_nextflow_config(mode)
    #script = "ete_build_dsl2.nf"
    cmd = [
        "nextflow", "run", script,
        "-ansi-log", "false",
	    "-profile", mode,  # Specify the profile based on the mode
        "--input", input_file,
        "--output", output_dir,
        "--aligner", aligner,
        "--trimmer", trimmer,
        "--tree_builder", tree_builder
    ]
    
    if resume:
        cmd.append("-resume")

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    return_code = process.poll()
    return return_code
    

def main():
    parser = argparse.ArgumentParser(description="Run Nextflow workflow.")
    parser.add_argument("--mode", default='local', choices=["local", "slurm"], required=True, help="Execution mode: local or slurm.")
    parser.add_argument("--partition", help="SLURM partition name (required if mode is slurm).")
    parser.add_argument("--time", default="1h", help="Time limit for SLURM jobs (only if mode is slurm).")
    parser.add_argument("--memory", default="4GB", help="Memory allocation for SLURM jobs (only if mode is slurm).")
    parser.add_argument("--cpus", type=int, default=4, help="Number of CPUs for SLURM jobs (only if mode is slurm).")
    parser.add_argument("--script", default="ete_build_dsl2.nf", help="Path to the Nextflow script to run.")
    parser.add_argument("--input", required=True, help="Input fasta file or directory.")
    parser.add_argument("--output", required=True, help="Output directory.")
    parser.add_argument("--aligner", default="mafft", help="Alignment tool.")
    parser.add_argument("--trimmer", default="trimal", help="Trimming tool.")
    parser.add_argument("--tree_builder", default="fasttree", help="Tree building tool.")
    parser.add_argument("--workflow", help="Select a predefined workflow.") #choices=list(PREDEFINED_WORKFLOWS.keys()),
    parser.add_argument("--resume", action="store_true", help="Resume from the last failed step.")

    args = parser.parse_args()
    
    # Validate SLURM-specific arguments
    if args.mode == "slurm" and not args.partition:
        parser.error("--partition is required when mode is slurm.")
    
    # If a predefined workflow is selected, override the tool choices
    if args.workflow:
        if PREDEFINED_WORKFLOWS.get(args.workflow):
            workflow_params = PREDEFINED_WORKFLOWS[args.workflow]
            args.aligner = workflow_params["aligner"]
            args.trimmer = workflow_params["trimmer"]
            args.tree_builder = workflow_params["tree_builder"]
        else:
            try:
                # parsing the workflow apps
                args.aligner, args.trimmer, args.tree_builder = args.workflow.split("-")
            except ValueError or IndexError:
                parser.error(f"Invalid workflow: {args.workflow}")

    run_nextflow(args.mode, args.input, args.output, args.aligner, args.trimmer, args.tree_builder, args.resume, args.script)

if __name__ == "__main__":
    main()
