#!/usr/bin/env python

import argparse
import os
import sys
import subprocess
import tempfile
import json

import src.mafft as mafft
import src.muscle as muscle
import src.clustalo as clustalo
import src.tcoffee as tcoffee
import src.famsa as famsa
import src.trimal as trimal
import src.trim_alg as trim_alg
import src.clipkit as clipkit
import src.fasttree as fasttree
import src.phyml as phyml


# Predefined workflows
PREDEFINED_WORKFLOWS = {
    "workflow1": {"aligner": "mafft", "trimmer": "trimal", "tree_builder": "fasttree"},
    "workflow2": {"aligner": "famsa", "trimmer": "trimal", "tree_builder": "fasttree"},
    "workflow3": {"aligner": "famsa", "trimmer": "trim_alg_v2", "tree_builder": "fasttree"},
    "ana-workflow": {"aligner": "hybrid", "trimmer": "trimal", "tree_builder": "fasttree"},
    # Add more predefined workflows as needed
}

ALIGNERS = ["mafft", "muscle", "tcoffee", "clustalo", "famsa"]
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


def remove_comment(line):
    """Remove comments from a line."""
    if '#' in line:
        line = line.split('#', 1)[0]
    return line.strip()


def convert_cfg_to_json(cfg_file, aligner, trimmer, tree_builder):
    """
    Convert a .cfg file to a JSON-like dictionary for the specified tools (aligner, trimmer, tree_builder).
    
    Parameters:
    - cfg_file (str): Path to the .cfg file.
    - aligner (str): The exact aligner to filter and include in the config.
    - trimmer (str): The exact trimmer to filter and include in the config.
    - tree_builder (str): The exact tree builder to filter and include in the config.
    
    Returns:
    - dict: A dictionary representing the filtered configuration for JSON.
    """
    config = {"aligner": {}, "trimmer": {}, "tree_builder": {}}
    current_section = None
    section_data = {}

    with open(cfg_file, 'r') as file:
        for line in file:
            line = remove_comment(line).strip()
            if not line or line.startswith("#"):
                continue

            if line.startswith("[") and line.endswith("]"):
                if current_section:
                    # Determine the correct parser based on the full section name
                    if current_section == aligner:
                        if section_data.get("_app") == "mafft":
                            mafft_config = mafft.parse_mafft_options(section_data)
                            config["aligner"]["mafft"] = mafft_config
                        elif section_data.get("_app") == "muscle":
                            config["aligner"]["muscle"] = muscle.parse_muscle_options(section_data)
                        elif section_data.get("_app") == "clustalo":
                            config["aligner"]["clustalo"] = clustalo.parse_clustalo_options(section_data)
                        elif section_data.get("_app") == "tcoffee":
                            config["aligner"]["tcoffee"] = tcoffee.parse_tcoffee_options(section_data)
                        elif section_data.get("_app") == "famsa":
                            config["aligner"]["famsa"] = famsa.parse_famsa_options(section_data)
                        else:
                            config["aligner"][section_data["_app"]] = section_data
                    elif current_section == trimmer:
                        if section_data.get("_app") == "trimal":
                            config["trimmer"]["trimal"] = trimal.parse_trimal_options(section_data)
                        elif section_data.get("_app") == "trim_alg_v2":
                            config["trimmer"]["trim_alg_v2"] = trim_alg.parse_trimalg_options(section_data)
                        elif section_data.get("_app") == "clipkit":
                            config["trimmer"]["clipkit"] = clipkit.parse_clipkit_options(section_data)
                        else:
                            config["trimmer"][section_data["_app"]] = section_data
                    elif current_section == tree_builder:
                        if section_data.get("_app") == "fasttree":
                            config["tree_builder"]["fasttree"] = fasttree.parse_fasttree_options(section_data)
                        elif section_data.get("_app") == "phyml":
                            config["tree_builder"]["phyml"] = phyml.parse_phyml_options(section_data)
                        else:
                            config["tree_builder"][section_data["_app"]] = section_data
                
                # Start a new section
                current_section = line[1:-1]  # Keep the full section name
                section_data = {}

            elif "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                value = remove_comment(value).strip()  # Remove comments from the value
                # Handle special case for empty strings
                if value == '""':
                    value = ""
                if value.isdigit():
                    value = int(value)
                elif value.replace('.', '', 1).isdigit():
                    value = float(value)
                elif value.lower() in ["true", "false"]:
                    value = value.lower() == "true"
                section_data[key] = value

        # Handle the last section
        if current_section:
            if current_section == aligner:
                if section_data.get("_app") == "mafft":
                    mafft_config = mafft.parse_mafft_options(section_data)
                    config["aligner"]["mafft"] = mafft_config
                elif section_data.get("_app") == "muscle":
                        config["aligner"]["muscle"] = muscle.parse_muscle_options(section_data)
                elif section_data.get("_app") == "clustalo":
                    config["aligner"]["clustalo"] = parse_clustalo_options(section_data)
                elif section_data.get("_app") == "tcoffee":
                            config["aligner"]["tcoffee"] = tcoffee.parse_tcoffee_options(section_data)
                elif section_data.get("_app") == "famsa":
                            config["aligner"]["famsa"] = famsa.parse_famsa_options(section_data)
                else:
                    config["aligner"][section_data["_app"]] = section_data
            elif current_section == trimmer:
                if section_data.get("_app") == "trimal":
                    config["trimmer"]["trimal"] = trimal.parse_trimal_options(section_data)
                elif section_data.get("_app") == "trim_alg_v2":
                    config["trimmer"]["trim_alg_v2"] = trim_alg.parse_trimalg_options(section_data)
                elif section_data.get("_app") == "clipkit":
                    config["trimmer"]["clipkit"] = clipkit.parse_clipkit_options(section_data)
                else:
                    config["trimmer"][section_data["_app"]] = section_data
            elif current_section == tree_builder:
                if section_data.get("_app") == "fasttree":
                    config["tree_builder"]["fasttree"] = fasttree.parse_fasttree_options(section_data)
                elif section_data.get("_app") == "phyml":
                    config["tree_builder"]["phyml"] = phyml.parse_phyml_options(section_data)
                else:
                    config["tree_builder"][section_data["_app"]] = section_data
    return config

def run_nextflow(mode, input_file, output_dir, aligner, trimmer, tree_builder, memory, threads, log_file, work_dir, workflow_config=None, resume=False, script="ete_build_dsl2.nf"):
    if workflow_config and workflow_config.endswith(".cfg"):
        cfg_json = convert_cfg_to_json(workflow_config, aligner, trimmer, tree_builder)
        json_file = workflow_config.replace(".cfg", ".json")
        
        with open(json_file, 'w') as out_json:
            json.dump(cfg_json, out_json, indent=4)
        workflow_config = json_file
    
    # Split aligner and tree_builder by underscore and take the first part
    aligner = aligner.split("_")[0]
    tree_builder = tree_builder.split("_")[0]
    
    # Handle the trimmer separately, preserving 'trim_alg_v2'
    if trimmer.startswith("trim_alg_v2"):
        trimmer = "trim_alg_v2"
    else:
        trimmer = trimmer.split("_")[0]
    
    cmd = [
        "nextflow", 
        "-C", "nextflow.config",
        "-log", log_file,
        "run", script,
        "--input", input_file,
        "--output", output_dir,
        "--aligner", aligner,
        "--trimmer", trimmer,
        "--tree_builder", tree_builder,
        "--memory", memory,
        "--thread", str(threads),
        "-work-dir", work_dir
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
    parser.add_argument("--trimmer", default="none", help="Trimming tool.")
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
