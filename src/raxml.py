def parse_raxml_options(section_data):
    """
    Parse RAxML options from the provided section data.
    
    Parameters:
    - section_data (dict): The section data extracted from the cfg file.
    
    Returns:
    - dict: A dictionary representing the RAxML configuration.
    """
    raxml_config = {
        "name": "raxml",
        "algorithm": section_data.get("algorithm", "d"),  # Default rapid hill-climbing algorithm
        "aa_model": section_data.get("aa_model", "PROTGAMMAJTT"),  # Default amino acid model
        "nt_model": section_data.get("nt_model", "GTRGAMMA"),  # Default nucleotide model
        "r_seed": section_data.get("r_seed", 31416),  # Default random seed
        "bootstrap": section_data.get("bootstrap", 0)  # Default to no bootstrap
    }

    return raxml_config