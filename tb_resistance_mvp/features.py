from .constants import HYDRO, VOLUME, CHARGE

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
