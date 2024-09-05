def parse_clipkit_options(section_data):
    """
    Parse Clipkit options from the provided section data.

    Parameters:
    - section_data (dict): The section data extracted from the cfg file.

    Returns:
    - dict: A dictionary representing the Clipkit configuration.
    """
    clipkit_config = {
        "name": "clipkit",
        "mode": section_data.get("mode", "smart-gap"),  # Default is smart-gap
        "gaps": section_data.get("gaps", 0.9),  # Default gap threshold is 0.9
        "codon": section_data.get("codon", False)  # Default codon mode is False
    }

    return clipkit_config
