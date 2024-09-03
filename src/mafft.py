def parse_mafft_options(section_data):
    """
    Parse MAFFT options from the provided section data.
    
    Parameters:
    - section_data (dict): The section data extracted from the cfg file.
    
    Returns:
    - dict: A dictionary representing the MAFFT configuration.
    """
    mafft_config = {
        "name": "mafft",
        #"mode": section_data.get("mode","auto"),
        "op": section_data.get("op", 1.53),
        "ep": section_data.get("ep", 0.123),
        "maxiterate": section_data.get("maxiterate", 0),
        "matrix": section_data.get("matrix", "").upper() if section_data.get("matrix") else "",
        "blosum_coefficient": section_data.get("blosum_coefficient", 62),
        "pam_coefficient": section_data.get("pam_coefficient", 80),
        "auto": section_data.get("auto", False),  # Automatically set to True if mode is "auto"
        "localpair": section_data.get("localpair", False),
        "globalpair": section_data.get("globalpair", False),
        "genafpair": section_data.get("genafpair", False),
        "retree": section_data.get("retree", None),
        "nofft": section_data.get("nofft", False),
        "parttree": section_data.get("parttree", False),
    }
    return mafft_config

if __name__ == "__main__":
    # Example usage
    section_example = {
        "_app": "mafft",
        #"mode": "linsi",
        "op": 1.53,
        "ep": 0.123,
        "maxiterate": 1000,
        "matrix": "",
        "blosum_coefficient": 62,
        "pam_coefficient": 80,
        "flag": "--localpair --maxiterate 1000"
    }
    config = parse_mafft_options(section_example)
    print(config)