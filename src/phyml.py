def parse_phyml_options(section_data):
    """
    Parse PhyML options from the provided section data.
    
    Parameters:
    - section_data (dict): The section data extracted from the cfg file.
    
    Returns:
    - dict: A dictionary representing the PhyML configuration.
    """
    phyml_config = {
        "name": "phyml",
        "aa_model": section_data.get("aa_model", "LG"),
        "nt_model": section_data.get("nt_model", "HKY85"),
        "pinv": section_data.get("pinv", "e"),  # Proportion of invariable sites
        "alpha": section_data.get("alpha", "e"),  # Gamma distribution parameter
        "nclasses": section_data.get("nclasses", 4),  # Number of rate categories
        "optimisation": section_data.get("optimisation", "tlr"),  # Optimisation parameters
        "frequencies": section_data.get("frequencies", "m"),  # Frequencies setting
        "bootstrap": section_data.get("bootstrap", -2),  # Bootstrap settings
        "tbe": section_data.get("tbe", False),  # TBE instead of FBP support
        "r_seed": section_data.get("r_seed", 123456),  # Random seed
    }
    
    return phyml_config

