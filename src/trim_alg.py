def parse_trimalg_options(section_data):
    """
    Parse trimalg options from the provided section data.
    
    Parameters:
    - section_data (dict): The section data extracted from the cfg file.
    
    Returns:
    - dict: A dictionary representing the Trimal configuration.
    """
    trimal_config = {
        "name": "trim_alg_v2",
        "min_res_abs": section_data.get("min_res_abs", None),  # Minimum number of residues (default: None)
        "min_res_percent": section_data.get("min_res_percent", None),  # Minimum percentage of residues (default: None)
    }
    return trimal_config