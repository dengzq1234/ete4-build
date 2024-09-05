def parse_fasttree_options(section_data):
    """
    Parse FastTree options from the provided section data.
    
    Parameters:
    - section_data (dict): The section data extracted from the .cfg file.
    
    Returns:
    - dict: A dictionary representing the FastTree configuration.
    """
    fasttree_config = {
        "name": "fasttree",
        "aa_model": section_data.get("aa_model", "JTT"),  # Amino acid model: LG, WAG, or JTT
        "nt_model": section_data.get("nt_model", "JC"),   # Nucleotide model: GTR or JC
        "gamma": section_data.get("gamma", False),        # Non-uniform evolutionary rates modeled by Gamma distribution
        "bootstrap": section_data.get("bootstrap", 1000), # Number of bootstrap replicates
        "pseudo": section_data.get("pseudo", False),      # Pseudo-likelihood support values
        "spr": section_data.get("spr", 4),                # Number of SPR rounds (minimum-evolution SPR moves)
        "mlacc": section_data.get("mlacc", 2),            # Rate categories for ML model of rate heterogeneity
        "slownni": section_data.get("slownni", False)     # Use slow NNI moves
    }
    return fasttree_config