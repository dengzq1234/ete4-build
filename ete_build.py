#!/usr/bin/env python

import argparse
import os
import sys
import subprocess
import tempfile
import json

# Predefined workflows
PREDEFINED_WORKFLOWS = {
    "workflow1": {"aligner": "mafft", "trimmer": "trimal", "tree_builder": "fasttree"},
    "workflow2": {"aligner": "famsa", "trimmer": "trimal", "tree_builder": "fasttree"},
    "workflow3": {"aligner": "famsa", "trimmer": "trim_alg_v2", "tree_builder": "fasttree"},
    "ana-workflow": {"aligner": "hybrid", "trimmer": "trimal", "tree_builder": "fasttree"},
    # Add more predefined workflows as needed
}

ALIGNERS = ["mafft", "muscle", "t_coffee", "clustalo", "famsa"]
TRIMMERS = ["trimal", "clipkit", "trim_alg_v2"]
TREE_BUILDERS = ["fasttree", "phyml", "raxml", "iqtree"]

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

def convert_cfg_to_json(cfg_file):
    """
    Convert a .cfg file to a JSON-like dictionary.
    """
    config = {"aligner": {}, "trimmer": {}, "tree_builder": {}}
    current_section = None
    section_data = {}

    with open(cfg_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if line.startswith("[") and line.endswith("]"):
                if current_section:
                    # Save the previous section to the correct category
                    app_type = section_data.get("_app")
                    if app_type in ALIGNERS:
                        config["aligner"][current_section] = section_data
                    elif app_type in TRIMMERS:
                        config["trimmer"][current_section] = section_data
                    elif app_type in TREE_BUILDERS:
                        config["tree_builder"][current_section] = section_data

                # Start a new section
                current_section = line[1:-1].replace("_default", "")
                section_data = {}

            elif "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip()
                if value.isdigit():
                    value = int(value)
                elif value.replace('.', '', 1).isdigit():
                    value = float(value)
                elif value.lower() in ["true", "false"]:
                    value = value.lower() == "true"
                section_data[key] = value

        # Save the last section
        if current_section:
            app_type = section_data.get("_app")
            if app_type in ALIGNERS:
                config["aligner"][current_section] = section_data
            elif app_type in TRIMMERS:
                config["trimmer"][current_section] = section_data
            elif app_type in TREE_BUILDERS:
                config["tree_builder"][current_section] = section_data

    return config

def run_nextflow(mode, input_file, output_dir, aligner, trimmer, tree_builder, memory, threads, log_file, work_dir, workflow_config=None, resume=False, script="ete_build_dsl2.nf"):
    # If a .cfg file is provided, convert it to a .json file

    if workflow_config and workflow_config.endswith(".cfg"):
        cfg_json = convert_cfg_to_json(workflow_config)
        json_file = workflow_config.replace(".cfg", ".json")
        with open(json_file, 'w') as out_json:
            json.dump(cfg_json, out_json, indent=4)
        workflow_config = json_file
    
    cmd = [
        "nextflow", 
        "-C", "nextflow.config",  # Specify the generated config file
        "-log", log_file,  # Add log file
        "run", script,
        "--input", input_file,
        "--output", output_dir,
        "--aligner", aligner,
        "--trimmer", trimmer,
        "--tree_builder", tree_builder,
        "--memory", memory,
        "--thread", str(threads),
        "-work-dir", work_dir  # Add work directory
    ]

    if workflow_config:
        cmd.extend(["--customConfig", workflow_config])
    
    if resume:
        cmd.append("-resume")

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
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
    parser.add_argument("--workflow", help="Select a predefined workflow.") #choices=list(PREDEFINED_WORKFLOWS.keys()),
    parser.add_argument("--resume", action="store_true", help="Resume from the last failed step.")
    parser.add_argument("--log", required=True, help="Log file location.")  # Log file argument
    parser.add_argument("--work-dir", required=True, help="Work directory location.")  # Work directory argument
    parser.add_argument("--config", help="Custom workflow config file.")

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

    run_nextflow(args.mode, args.input, args.output, args.aligner, args.trimmer, args.tree_builder, args.memory, args.cpus, args.log, args.work_dir, args.config, args.resume, args.script)

if __name__ == "__main__":
    main()
