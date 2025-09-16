AAS = list("ACDEFGHIKLMNPQRSTVWY")

# Kyte-Doolittle hydropathy index (approx)
HYDRO = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
    'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
    'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

# Approximate side-chain volumes (Å^3), small set to create deltas
VOLUME = {
    'G': 48, 'A': 67, 'S': 73, 'P': 90, 'V': 105, 'T': 93, 'C': 86,
    'I': 124, 'L': 124, 'N': 96, 'D': 91, 'Q': 114, 'K': 135, 'E': 109,
    'M': 124, 'H': 118, 'F': 135, 'R': 148, 'Y': 141, 'W': 163
}

# Charge groups at physiological pH (very rough)
CHARGE = {
    'D': -1, 'E': -1, 'K': +1, 'R': +1, 'H': +0.1  # H partially positive
}

# Drug panels by gene hint (you can extend this)
GENE_DRUGS = {
    "rpoB": ["RIF", "RFB"],
    "katG": ["INH"],
    "pncA": ["PZA"],
    "embB": ["EMB"],
    "gyrA": ["FQ"],
    "gyrB": ["FQ"],
    "rrs": ["AMK", "KAN", "CAP"],
    "eis": ["KAN"],
    "rplC": ["LZD"],
    "rplD": ["LZD"],
    "atpE": ["BDQ"],
    "ddn": ["DLM", "PMD"],
}
DEFAULT_DRUGS = ["RIF", "INH", "PZA", "EMB", "FQ", "AMK", "KAN", "CAP", "LZD", "BDQ", "DLM", "PMD"]
