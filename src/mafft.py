# mafft.py

def parse_mafft_options(aligner_option):
    """
    Parse the MAFFT aligner option and return the appropriate JSON configuration.
    
    Parameters:
    - aligner_option (str): The MAFFT option passed via command-line (e.g., "mafft_default", "mafft_linsi")
    
    Returns:
    - dict: A dictionary representing the MAFFT configuration to be used in the JSON config.
    """
    # Base configuration for MAFFT
    mafft_base_config = {
        "name": "mafft",
        "mode": "auto",
        "op": 1.53,
        "ep": 0.123,
        "maxiterate": 0,
        "matrix": "",
        "blosum_coefficient": 62,
        "pam_coefficient": 80,
        "methods": {}
    }

    # Method-specific configurations
    method_specific_configs = {
        "mafft_default": {
            "methods": {
                "auto": {"flag": "--auto"}
            }
        },
        "mafft_linsi": {
            "methods": {
                "linsi": {"flag": "--localpair --maxiterate 1000"}
            }
        },
        "mafft_ginsi": {
            "methods": {
                "ginsi": {"flag": "--globalpair --maxiterate 1000"}
            }
        },
        "mafft_einsi": {
            "methods": {
                "einsi": {"flag": "--ep 0 --genafpair --maxiterate 1000"}
            }
        },
        "mafft_fftnsi": {
            "methods": {
                "fftnsi": {"flag": "--retree 2 --maxiterate 2"}
            }
        },
        "mafft_fftnsi_max": {
            "methods": {
                "fftnsi_max": {"flag": "--retree 2 --maxiterate 1000"}
            }
        },
        "mafft_fftns": {
            "methods": {
                "fftns": {"flag": "--retree 2 --maxiterate 0"}
            }
        },
        "mafft_fftns1": {
            "methods": {
                "fftns1": {"flag": "--retree 1 --maxiterate 0"}
            }
        },
        "mafft_nwnsi": {
            "methods": {
                "nwnsi": {"flag": "--retree 2 --maxiterate 2 --nofft"}
            }
        },
        "mafft_nwns": {
            "methods": {
                "nwns": {"flag": "--retree 2 --maxiterate 0 --nofft"}
            }
        },
        "mafft_nwns_parttree": {
            "methods": {
                "nwns_parttree": {"flag": "--retree 1 --maxiterate 0 --nofft --parttree"}
            }
        }
    }

    # Apply the specific method configuration to the base config
    if aligner_option in method_specific_configs:
        for key, value in method_specific_configs[aligner_option].items():
            if key == "methods":
                mafft_base_config["methods"].update(value)
            else:
                mafft_base_config[key] = value
    else:
        raise ValueError(f"Unknown MAFFT option: {aligner_option}")

    return mafft_base_config


if __name__ == "__main__":
    # Example usage
    import sys
    aligner_option = sys.argv[1]  # e.g., "mafft_linsi"
    mafft_config = parse_mafft_options(aligner_option)
    print(mafft_config)
