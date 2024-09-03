def parse_muscle_options(section_data):
    """
    Parse MUSCLE options from the provided section data.
    
    Parameters:
    - section_data (dict): The section data extracted from the cfg file.
    
    Returns:
    - dict: A dictionary representing the MUSCLE configuration.
    """
    muscle_config = {
        "name": "muscle",
        #"mode": section_data.get("mode", "align"),
        "perturb": section_data.get("perturb", 0),
        "perm": section_data.get("perm", "").lower(),  # Keeping perm as lower-case since it's more natural for a command-line flag
        "replicates": section_data.get("replicates", None),  # Optional, specific to certain modes
        "stratified": section_data.get("stratified", False),
        "diversified": section_data.get("diversified", False),
        "consiters": section_data.get("consiters", 2),  # Optional, defaulting to 2
        "refineiters": section_data.get("refineiters", 100),  # Optional, defaulting to 100
    }

    return muscle_config
