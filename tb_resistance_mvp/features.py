from .constants import HYDRO, VOLUME, CHARGE

# Codon table for translation
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def translate_codon(codon: str) -> str:
    """Translate a single codon to amino acid"""
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3:
        return 'X'
    return CODON_TABLE.get(codon, 'X')

def physchem_features(wt: str, mut: str, pos: int, seq_len: int):
    # Hydropathy change
    dh = (HYDRO.get(mut, 0.0) - HYDRO.get(wt, 0.0))
    # Volume change (normalized)
    dv = (VOLUME.get(mut, 110) - VOLUME.get(wt, 110)) / 110.0
    # Charge change
    dc = (CHARGE.get(mut, 0.0) - CHARGE.get(wt, 0.0))
    # Positional features
    rel_pos = pos / seq_len
    center_dist = abs(rel_pos - 0.5) * 2  # 0 center, 1 at ends
    # Simple motif penalty: mutations to Proline or Glycine in middle region are more disruptive
    motif_pen = 1.0 if (mut in ("P","G") and 0.2 < rel_pos < 0.8) else 0.0
    return {
        "d_hydro": dh,
        "d_vol": dv,
        "d_charge": dc,
        "rel_pos": rel_pos,
        "center_dist": center_dist,
        "motif_pen": motif_pen
    }
