def parse_clustalo_options(section_data):
    """
    Parse Clustal Omega (ClustalO) options from the provided section data.
    
    Parameters:
    - section_data (dict): The section data extracted from the cfg file.
    
    Returns:
    - dict: A dictionary representing the ClustalO configuration.
    """
    clustalo_config = {
        "name": "clustalo",
        "dealign": section_data.get("dealign", False),
        "full": section_data.get("full", False),
        "full_iter": section_data.get("full_iter", False),
        "iterations": section_data.get("iterations", None),
        "max_guidetree_iterations": section_data.get("max_guidetree_iterations", None),
        "max_hmm_iterations": section_data.get("max_hmm_iterations", None)
    }

    return clustalo_config