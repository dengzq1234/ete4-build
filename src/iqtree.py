def parse_iqtree_options(section_data):
    """
    Parse IQ-TREE options from the provided section data.
    
    Parameters:
    - section_data (dict): The section data extracted from the cfg file.
    
    Returns:
    - dict: A dictionary representing the IQ-TREE configuration.
    """
    iqtree_config = {
        "name": "iqtree",
        "alrt": section_data.get("alrt", 0),  # Approximate likelihood ratio test support (default: 0)
        "ufboot": section_data.get("ufboot", None),  # Ultrafast bootstrap replicates (optional)
        "seed": section_data.get("seed", 31416),  # Default seed
        "model": section_data.get("model", "TEST"),  # Default to TEST
        "tbe": section_data.get("tbe", False),  # Ultrafast bootstrap with transfer bootstrap expectation
        "st": section_data.get("st", None)  # Codon substitution models (optional)
    }

    return iqtree_config