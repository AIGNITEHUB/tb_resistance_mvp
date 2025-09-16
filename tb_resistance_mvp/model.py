from math import tanh

# Very simple, deterministic "scorer" per drug using feature weights
# This is a placeholder for a real ML model; output is 0..1

DRUG_WEIGHTS = {
    # drug: (w_dh, w_dv, w_dc, w_center, w_motif, bias)
    "RIF": ( 0.6,  0.4,  0.2,  0.1, 0.5,  0.0),
    "RFB": ( 0.6,  0.4,  0.2,  0.1, 0.5, -0.05),
    "INH": ( 0.2,  0.3,  0.6,  0.1, 0.4,  0.0),
    "PZA": ( 0.3,  0.4,  0.3,  0.2, 0.3, -0.05),
    "EMB": ( 0.2,  0.5,  0.2,  0.2, 0.2, -0.05),
    "FQ":  ( 0.4,  0.2,  0.3,  0.3, 0.4,  0.0),
    "AMK": ( 0.3,  0.2,  0.5,  0.2, 0.3,  0.0),
    "KAN": ( 0.3,  0.2,  0.6,  0.2, 0.3,  0.0),
    "CAP": ( 0.2,  0.2,  0.5,  0.2, 0.3, -0.05),
    "LZD": ( 0.2,  0.2,  0.3,  0.3, 0.4, -0.05),
    "BDQ": ( 0.3,  0.5,  0.2,  0.2, 0.3,  0.0),
    "DLM": ( 0.3,  0.4,  0.2,  0.3, 0.3, -0.05),
    "PMD": ( 0.3,  0.4,  0.2,  0.3, 0.3, -0.05),
}

def sigmoid(x: float) -> float:
    return 1 / (1 + pow(2.718281828, -x))

def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))

class HeuristicModel:
    def __init__(self, drugs):
        self.drugs = drugs

    def predict(self, feats: dict) -> dict:
        out = {}
        for d in self.drugs:
            w = DRUG_WEIGHTS.get(d, (0.3,0.3,0.3,0.2,0.2,0.0))
            score = (
                w[0]*feats["d_hydro"] +
                w[1]*feats["d_vol"]   +
                w[2]*feats["d_charge"]+
                w[3]*feats["center_dist"] +
                w[4]*feats["motif_pen"] +
                w[5]
            )
            # squash to 0..1; tanh to smooth extremes
            p = (tanh(score) + 1)/2
            out[d] = clamp01(p)
        return out
