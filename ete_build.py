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
    "workflow3": {"aligner": "famsa", "trimmer": "trim_alg_v2", "tree_builder": "fasttree"},
    "ana-workflow": {"aligner": "hybrid", "trimmer": "trimal", "tree_builder": "fasttree"},
    # Add more predefined workflows as needed
}

# def generate_nextflow_config(execution_mode, partition=None, time=None, memory=None, cpus=None):
#     if execution_mode == "slurm":
#         config_content = f"""
#             profiles {{
#                 slurm {{
#                     process {{
#                         executor = 'slurm'
#                         queue = '{partition}'
#                         time = '{time}'
#                         memory = '{memory}'
#                         cpus = {cpus}
#                     }}
#                 }}
#             }}
#             """
#     elif execution_mode == "local":
#         config_content = """
# profiles {
#     local {
#         process {
#             executor = 'local'
#         }
#     }
# }
# """
#     else:
#         raise ValueError(f"Unsupported execution mode: {execution_mode}")

#     with open("nextflow.config", "w") as f:
#         f.write(config_content)

def generate_nextflow_config(args):
    """
    Generate the nextflow.config file based on the provided arguments.
    """
    # Base configuration content
    config_content = f"""
params.input = '{args.input}'
params.output = '{args.output}'
params.thread = {args.cpus}
params.aligner = '{args.aligner}'
params.trimmer = '{args.trimmer}'
params.tree_builder = '{args.tree_builder}'
params.memory = '{args.memory}'
params.time = '{args.time}'
"""

    # Add profile-specific configurations
    if args.mode == "slurm":
        config_content += f"""
profiles {{
    slurm {{
        process {{
            executor = 'slurm'
            queue = '{args.partition}'
            time = '{args.time}'
            memory = '{args.memory}'
            cpus = {args.cpus}
        }}
    }}
}}
"""
    elif args.mode == "local":
        config_content += """
profiles {
    local {
        process {
            executor = 'local'
        }
    }
}
"""
    else:
        raise ValueError(f"Unsupported execution mode: {args.mode}")

    # Write to the nextflow.config file
    with open("nextflow.config", "w") as f:
        f.write(config_content)

def run_nextflow(mode, input_file, output_dir, aligner, trimmer, tree_builder, memory, threads, resume=False, script="ete_build_dsl2.nf"):
    #script = "ete_build_dsl2.nf"
    # cmd = [
    #     "nextflow", "run", script,
    #     "-ansi-log", "false",
	#     "-profile", mode,  # Specify the profile based on the mode
    #     "--input", input_file,
    #     "--output", output_dir,
    #     "--aligner", aligner,
    #     "--trimmer", trimmer,
    #     "--tree_builder", tree_builder,
    #     "--memory", memory,
    #     "--thread", str(threads),
    # ]
    
    cmd = [
        "nextflow", 
        "-C", "nextflow.config",  # Specify the generated config file
        "run", script,
       
        
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

    # Check for errors
    if return_code != 0:
        print(f"Error: Nextflow script {script} exited with code {return_code}", file=sys.stderr)

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
    parser.add_argument("--workflow", choices=list(PREDEFINED_WORKFLOWS.keys()), help="Select a predefined workflow.")
    parser.add_argument("--resume", action="store_true", help="Resume from the last failed step.")

    args = parser.parse_args()
    
    # Validate SLURM-specific arguments
    if args.mode == "slurm" and not args.partition:
        parser.error("--partition is required when mode is slurm.")
    
    # If a predefined workflow is selected, override the tool choices
    if args.workflow:
        workflow_params = PREDEFINED_WORKFLOWS[args.workflow]
        args.aligner = workflow_params["aligner"]
        args.trimmer = workflow_params["trimmer"]
        args.tree_builder = workflow_params["tree_builder"]

    # Generate the Nextflow config AFTER setting the workflow-specific parameters
    generate_nextflow_config(args)

    run_nextflow(args.mode, args.input, args.output, args.aligner, args.trimmer, args.tree_builder, args.memory, args.cpus, args.resume, args.script)

if __name__ == "__main__":
    main()
