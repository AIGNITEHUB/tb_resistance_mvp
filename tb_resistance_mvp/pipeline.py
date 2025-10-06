from typing import List, Dict
from .constants import AAS, GENE_DRUGS, DEFAULT_DRUGS
from .features import physchem_features
from .model import HeuristicModel
from statistics import mean

# Import demo version (no RDKit needed!)
try:
    from .drug_filtering_demo import MutationGroupClassifier, FakeMoleculeDatabase, FilterCriteriaHelper
    DEMO_MODE = True
except ImportError:
    # Fallback if demo module not available
    DEMO_MODE = False


def _window_scores(heatmap_drug: dict, start: int, end: int, w: int = 12):
    """Average max P(resistance) per sliding window (lower = better for design)."""
    scores = []
    for s in range(start, end - w + 2):
        e = s + w - 1
        vals = [heatmap_drug[p] for p in range(s, e + 1)]
        scores.append({"start": s, "end": e, "w": w, "avg_escape": float(mean(vals))})
    return scores


def _suggest_scaffolds(gene_hint: str, drugs: list, region: dict):
    """Heuristic scaffold suggestions for in-silico screening."""
    base = [
        {"name": "bicyclic heteroaromatic", "smiles": "c1ncccc1c2ncccc2", "note": "rigid, N-rich"},
        {"name": "quinazolinone-like", "smiles": "O=C1N=CNc2ccccc12", "note": "dual H-bond capability"},
        {"name": "benzoxazole-like", "smiles": "O1c2ccccn2C=N1", "note": "small, moderate hydrophobicity"},
    ]
    
    if gene_hint.lower() in ("rpob",):
        props = {"target_drugs": ["RIF", "RFB"], "mw": "350–520", "cLogP": "2.0–4.0", "hbd/hba": "≤1 / 3–6"}
    elif gene_hint.lower() in ("pnca",):
        props = {"target_drugs": ["PZA"], "mw": "200–380", "cLogP": "1.0–3.0", "hbd/hba": "0–1 / 2–5"}
    elif gene_hint.lower() in ("katg", "inha"):
        props = {"target_drugs": ["INH"], "mw": "250–420", "cLogP": "1.5–3.5", "hbd/hba": "1–2 / 3–6"}
    else:
        props = {"target_drugs": drugs, "mw": "250–500", "cLogP": "1.5–4.0", "hbd/hba": "≤2 / 3–6"}
    
    return {"props": props, "seed_scaffolds": base, "region_hint": region}


def _filter_high_risk_mutations(records: list, cutoff: float = 0.7) -> list:
    """Filter mutations with high resistance probability."""
    high_risk = []
    for r in records:
        max_prob = max(r["p_resist"].values())
        if max_prob >= cutoff:
            high_risk.append({
                **r,
                "max_prob": max_prob,
                "risk_drugs": [d for d, p in r["p_resist"].items() if p >= cutoff]
            })
    return high_risk


def _group_mutations_by_class(high_risk_mutations: list) -> dict:
    """Group high-risk mutations by physicochemical class."""
    grouped = {}
    for mut in high_risk_mutations:
        group = mut["group"]
        if group not in grouped:
            grouped[group] = []
        grouped[group].append(mut)
    return grouped


def infer_drugs(gene_hint: str) -> List[str]:
    if not gene_hint:
        return DEFAULT_DRUGS
    return GENE_DRUGS.get(gene_hint, DEFAULT_DRUGS)


def mutational_scan(
    sequence: str,
    region_start: int,
    region_end: int,
    gene_hint: str = "",
    enable_drug_filtering: bool = True,
    resistance_cutoff: float = 0.7,
    demo_mode: bool = True  # NEW: Enable demo mode by default
) -> Dict:
    """
    Perform mutational scanning with optional drug candidate filtering
    
    Args:
        sequence: Protein sequence
        region_start: Start position (1-indexed)
        region_end: End position (1-indexed)
        gene_hint: Gene name hint (rpoB, katG, etc.)
        enable_drug_filtering: Enable 6-group classification and filtering
        resistance_cutoff: Probability threshold for high-risk classification
        demo_mode: Use fake molecules for demo (no RDKit needed)
    
    Returns:
        Dictionary with mutation matrix, heatmap, and drug suggestions
    """
    sequence = sequence.replace("\n", "").replace(" ", "").strip().upper()
    L = len(sequence)
    if L == 0:
        raise ValueError("Empty sequence")
    if not (1 <= region_start <= L and 1 <= region_end <= L and region_start <= region_end):
        raise ValueError("Invalid region")
    
    drugs = infer_drugs(gene_hint)
    model = HeuristicModel(drugs)
    records = []
    
    for pos1 in range(region_start, region_end + 1):
        wt = sequence[pos1 - 1]
        for mut in AAS:
            if mut == wt:
                continue
            
            feats = physchem_features(wt, mut, pos1, L)
            probs = model.predict(feats)
            
            # Simple mechanism tag (legacy)
            mech = "binding_site" if feats["d_hydro"] > 1.5 else (
                "stability" if abs(feats["d_vol"]) > 0.4 else (
                    "charge_effect" if abs(feats["d_charge"]) > 0 else "neutralish"
                )
            )
            
            # NEW: Classify into 6 groups
            group = "General"
            if enable_drug_filtering and DEMO_MODE:
                group = MutationGroupClassifier.classify(feats, wt, mut)
            
            records.append({
                "pos": pos1,
                "wt": wt,
                "mut": mut,
                "features": feats,
                "p_resist": probs,
                "mechanism": mech,
                "group": group
            })
    
    # Aggregate heatmap
    heatmap = {d: {} for d in drugs}
    for d in drugs:
        for pos1 in range(region_start, region_end + 1):
            maxp = max([r["p_resist"][d] for r in records if r["pos"] == pos1], default=0.0)
            heatmap[d][pos1] = maxp
    
    # Select low-escape regions
    main_drug = drugs[0]
    windows = _window_scores(heatmap[main_drug], region_start, region_end, w=12)
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
            "rationale": "low escape probability"
        })
        used.append(w)
        if len(design_regions) >= 3:
            break
    
    # Scaffold suggestions
    design = _suggest_scaffolds(gene_hint, drugs, design_regions[0] if design_regions else {
        "start": region_start, "end": region_end})
    
    result = {
        "matrix": records,
        "heatmap": heatmap,
        "drugs": drugs,
        "region": {"start": region_start, "end": region_end},
        "length": L,
        "design_regions": design_regions,
        "design": design
    }
    
    # NEW: Drug filtering analysis
    if enable_drug_filtering and DEMO_MODE:
        high_risk = _filter_high_risk_mutations(records, resistance_cutoff)
        grouped = _group_mutations_by_class(high_risk)
        
        # Generate filter summary and example molecules
        filter_summary = {}
        example_molecules = {}
        
        for group in grouped.keys():
            filter_summary[group] = {
                "name": FilterCriteriaHelper.CRITERIA.get(group, {}).get("name", group),
                "description": FilterCriteriaHelper.get_description(group),
                "mutation_count": len(grouped[group]),
                "example_mutations": [
                    f"{m['wt']}{m['pos']}{m['mut']}"
                    for m in grouped[group][:5]
                ]
            }
            
            # NEW: Generate example molecules for each group (DEMO)
            if demo_mode:
                molecules = FakeMoleculeDatabase.search(group, max_results=10)
                example_molecules[group] = molecules
        
        result["drug_filtering"] = {
            "enabled": True,
            "demo_mode": demo_mode,
            "cutoff": resistance_cutoff,
            "high_risk_count": len(high_risk),
            "groups": grouped,
            "filter_summary": filter_summary,
            "example_molecules": example_molecules  # NEW
        }
    
    return result