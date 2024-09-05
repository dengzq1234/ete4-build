def parse_famsa_options(section_data):
    """
    Parse FAMSA options from the provided section data.
    
    Parameters:
    - section_data (dict): The section data extracted from the cfg file.
    
    Returns:
    - dict: A dictionary representing the FAMSA configuration.
    """
    famsa_config = {
        "name": "famsa",
        "gt": section_data.get("gt", None),  # Guide tree: "", sl, upgma, or nj, or import a tree in newick format
        "medoidtree": section_data.get("medoidtree", False),  # Use medoid guide tree, defaults to False
        "refine_mode": section_data.get("refine_mode", "auto"),  # auto, off, or on (default: auto)
        "r": section_data.get("r", 100),  # Refinement iterations (default is 0 if not provided)
        "go": section_data.get("go", None),  # Gap open penalty
        "ge": section_data.get("ge", None),  # Gap extension penalty
        "tgo": section_data.get("tgo", None),  # Terminal gap open penalty
        "tge": section_data.get("tge", None),  # Terminal gap extension penalty
        "gsd": section_data.get("gsd", None),  # Gap cost scaler div-term
        "gsl": section_data.get("gsl", None),  # Gap cost scaler log-term
        "dgr": section_data.get("dgr", False),  # Disable gap cost rescaling, defaults to False
        "dgo": section_data.get("dgo", False),  # Disable gap open optimization, defaults to False
        "dsp": section_data.get("dsp", False)  # Disable sum of pairs optimization, defaults to False
    }
    
    return famsa_config
