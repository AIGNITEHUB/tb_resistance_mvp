from typing import List, Dict
from .constants import AAS, GENE_DRUGS, DEFAULT_DRUGS
from .features import physchem_features
from .model import HeuristicModel
from statistics import mean

def _window_scores(heatmap_drug: dict, start: int, end: int, w: int = 12):
    """Average max P(resistance) per sliding window (lower = tốt cho thiết kế)."""
    scores = []
    for s in range(start, end - w + 2):
        e = s + w - 1
        vals = [heatmap_drug[p] for p in range(s, e + 1)]
        scores.append({"start": s, "end": e, "w": w, "avg_escape": float(mean(vals))})
    return scores

def _suggest_scaffolds(gene_hint: str, drugs: list, region: dict):
    """Heuristic gợi ý scaffold + thuộc tính đích cho sàng lọc in-silico (không RDKit)."""
    # gợi ý mặc định
    base = [
        {"name": "bicyclic heteroaromatic", "smiles": "c1ncccc1c2ncccc2", "note":"khung cứng, giàu N"},
        {"name": "quinazolinone-like",     "smiles": "O=C1N=CNc2ccccc12",   "note":"khả năng H-bond đôi"},
        {"name": "benzoxazole-like",       "smiles": "O1c2ccccn2C=N1",      "note":"nhỏ, kỵ nước vừa"},
    ]
    # tuỳ gene → điều chỉnh độ phân cực kích thước
    if gene_hint.lower() in ("rpoB",):
        props = {"target_drugs": ["RIF","RFB"], "mw": "350–520", "cLogP": "2.0–4.0", "hbd/hba": "≤1 / 3–6"}
    elif gene_hint.lower() in ("pncA",):
        props = {"target_drugs": ["PZA"], "mw": "200–380", "cLogP": "1.0–3.0", "hbd/hba": "0–1 / 2–5"}
    elif gene_hint.lower() in ("katG","inhA"):
        props = {"target_drugs": ["INH"], "mw": "250–420", "cLogP": "1.5–3.5", "hbd/hba": "1–2 / 3–6"}
    else:
        props = {"target_drugs": drugs, "mw": "250–500", "cLogP": "1.5–4.0", "hbd/hba": "≤2 / 3–6"}
    return {"props": props, "seed_scaffolds": base, "region_hint": region}

def infer_drugs(gene_hint: str) -> List[str]:
    if not gene_hint:
        return DEFAULT_DRUGS
    return GENE_DRUGS.get(gene_hint, DEFAULT_DRUGS)

def mutational_scan(sequence: str, region_start: int, region_end: int, gene_hint: str = "") -> Dict:
    sequence = sequence.replace("\n", "").replace(" ", "").strip().upper()
    L = len(sequence)
    if L == 0:
        raise ValueError("Empty sequence")
    if not (1 <= region_start <= L and 1 <= region_end <= L and region_start <= region_end):
        raise ValueError("Invalid region")
    drugs = infer_drugs(gene_hint)
    model = HeuristicModel(drugs)
    records = []
    for pos1 in range(region_start, region_end+1):  # 1-indexed inclusive
        wt = sequence[pos1-1]
        for mut in AAS:
            if mut == wt:
                continue
            feats = physchem_features(wt, mut, pos1, L)
            probs = model.predict(feats)
            # simple mechanism tag
            mech = "binding_site" if feats["d_hydro"] > 1.5 else ("stability" if abs(feats["d_vol"])>0.4 else "charge_effect" if abs(feats["d_charge"])>0 else "neutralish")
            records.append({
                "pos": pos1,
                "wt": wt,
                "mut": mut,
                "features": feats,
                "p_resist": probs,
                "mechanism": mech
            })
    # aggregate a simple heatmap per drug: max prob per position
    heatmap = {d: {} for d in drugs}
    for d in drugs:
        for pos1 in range(region_start, region_end+1):
            maxp = max([r["p_resist"][d] for r in records if r["pos"]==pos1], default=0.0)
            heatmap[d][pos1] = maxp

    # --- NEW: chọn vùng "ít escape" ưu tiên (ví dụ dùng thuốc đầu trong panel) ---
    main_drug = drugs[0]
    windows = _window_scores(heatmap[main_drug], region_start, region_end, w=12)
    # xếp theo avg_escape tăng dần, lọc 3 vùng đầu và loại trùng chồng lấn
    windows.sort(key=lambda x: x["avg_escape"])
    design_regions = []
    used = []
    for w in windows:
        if any(not (w["end"] < u["start"] or w["start"] > u["end"]) for u in used):
            continue
        design_regions.append({
            "start": w["start"], "end": w["end"],
            "avg_escape": round(w["avg_escape"], 4),
            "drug": main_drug,
            "rationale": "escape thấp (trung bình)"
        })
        used.append(w)
        if len(design_regions) >= 3:
            break

    # --- NEW: đề xuất scaffold/thuộc tính ---
    design = _suggest_scaffolds(gene_hint, drugs, design_regions[0] if design_regions else {"start":region_start,"end":region_end})

    return {
        "matrix": records,
        "heatmap": heatmap,
        "drugs": drugs,
        "region": {"start": region_start, "end": region_end},
        "length": L,
        "design_regions": design_regions,
        "design": design
    }
