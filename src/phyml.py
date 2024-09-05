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
        "datatype": section_data.get("datatype", "aa"),
    }
    return phyml_config

