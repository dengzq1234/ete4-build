def parse_mrbayes_options(section_data):
    """
    Parse MrBayes options from the provided section data.

    Parameters:
    - section_data (dict): The section data extracted from the cfg file.

    Returns:
    - dict: A dictionary representing the MrBayes configuration.
    """
    mrbayes_config = {
        "name": "mrbayes",
        "ngen": int(section_data.get("ngen", 100000)),  # Number of generations (default: 100000)
        "nchains": int(section_data.get("nchains", 4)),  # Number of chains (default: 4)
        "nruns": int(section_data.get("nruns", 1)),  # Number of runs (default: 1)
        "nst": section_data.get("nst", 6),  # Substitution model for DNA (default: 6)
        "rates": section_data.get("rates", "invgamma"),  # Rates variation for DNA (default: invgamma)
        "aamodelpr": section_data.get("aamodelpr", "fixed(wag)"),  # Amino acid model prior (default: fixed(wag))
        "diagnfreq": int(section_data.get("diagnfreq", 5000)),  # Frequency of diagnostics (default: 5000)
        "samplefreq": int(section_data.get("samplefreq", 500)),  # Frequency of sampling (default: 500)
        "printfreq": int(section_data.get("printfreq", 1000)),  # Frequency of printing (default: 1000)
        "burninfrac": float(section_data.get("burninfrac", 0.25)),  # Fraction of samples to discard (default: 0.25)
        "append": "yes" if section_data.get("append", "False").lower() == "true" else "no",  # Append to the previous analysis (default: no)
        "stoprule": "yes" if section_data.get("stoprule", "False").lower() == "true" else "no",  # Use stop rule (default: no)
        "seed": int(section_data.get("seed", 123456)),  # Random seed (default: 123456)
        "swapseed": int(section_data.get("swapseed", 123456)),  # Swap seed (default: 123456)
    }

    return mrbayes_config
