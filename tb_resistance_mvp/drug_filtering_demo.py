"""
Demo version of drug filtering - NO RDKit required
Uses simple heuristics and fake molecule database
"""

from typing import Dict, List
import random


class MutationGroupClassifier:
    """Classify mutations into 6 groups based on physicochemical features"""
    
    THRESHOLDS = {
        "d_hydro_h": -2.5,
        "d_hydro_p": 3.0,
        "d_charge": 1.0,
        "d_vol": 0.6,
        "motif_pen": 0.7,
        "rel_pos": 0.8
    }
    
    @classmethod
    def classify(cls, features: dict, wt: str, mut: str) -> str:
        """Classify mutation into one of 6 groups"""
        dh = features["d_hydro"]
        dv = features["d_vol"]
        dc = features["d_charge"]
        mp = features["motif_pen"]
        rp = features["rel_pos"]
        
        # Priority order
        if mp >= cls.THRESHOLDS["motif_pen"] or rp >= cls.THRESHOLDS["rel_pos"]:
            return "M"
        
        if abs(dc) >= cls.THRESHOLDS["d_charge"]:
            return "C+" if dc > 0 else "C-"
        
        aromatic_aas = {'F', 'Y', 'W', 'H'}
        if wt in aromatic_aas and mut not in aromatic_aas:
            return "A"
        
        if abs(dv) >= cls.THRESHOLDS["d_vol"]:
            return "V"
        
        if dh <= cls.THRESHOLDS["d_hydro_h"]:
            return "H"
        
        if dh >= cls.THRESHOLDS["d_hydro_p"]:
            return "P"
        
        return "General"


class FakeMoleculeDatabase:
    """
    Fake molecule database for demo purposes
    In real version, this would be replaced by ZINC + RDKit
    """
    
    # Fake SMILES strings (simplified representations)
    FAKE_MOLECULES = {
        "H": [  # Hydrophobic
            {"name": "Naphthalene-like", "smiles": "c1ccc2ccccc2c1", "mw": 128, "logp": 3.4, "tpsa": 20},
            {"name": "Biphenyl-like", "smiles": "c1ccc(cc1)c2ccccc2", "mw": 154, "logp": 4.1, "tpsa": 15},
            {"name": "Indole-analog", "smiles": "c1ccc2c(c1)[nH]cc2", "mw": 117, "logp": 2.9, "tpsa": 28},
        ],
        "P": [  # Polar
            {"name": "Glucose-like", "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O", "mw": 180, "logp": -2.1, "tpsa": 110},
            {"name": "Amino-acid-analog", "smiles": "NC(C(=O)O)CO", "mw": 105, "logp": -3.2, "tpsa": 83},
            {"name": "Glycerol-derivative", "smiles": "OCC(O)CO", "mw": 92, "logp": -1.8, "tpsa": 60},
        ],
        "C+": [  # Need anionic (for positive charge site)
            {"name": "Benzoic-acid", "smiles": "c1ccc(cc1)C(=O)[O-]", "mw": 122, "logp": 1.9, "tpsa": 40, "charge": -1},
            {"name": "Sulfonic-acid", "smiles": "c1ccccc1S(=O)(=O)[O-]", "mw": 158, "logp": 0.8, "tpsa": 60, "charge": -1},
            {"name": "Phosphate-ester", "smiles": "CCOP(=O)([O-])[O-]", "mw": 140, "logp": -1.2, "tpsa": 75, "charge": -2},
        ],
        "C-": [  # Need cationic (for negative charge site)
            {"name": "Phenethylamine", "smiles": "c1ccc(cc1)CC[NH3+]", "mw": 122, "logp": 1.4, "tpsa": 26, "charge": +1},
            {"name": "Lysine-analog", "smiles": "NCCCC[NH3+]", "mw": 103, "logp": -2.8, "tpsa": 52, "charge": +1},
            {"name": "Guanidinium", "smiles": "NC(=[NH2+])N", "mw": 59, "logp": -1.5, "tpsa": 88, "charge": +1},
        ],
        "V": [  # Volume change
            {"name": "Ethanol-small", "smiles": "CCO", "mw": 46, "logp": -0.3, "tpsa": 20},
            {"name": "Adamantane-bulky", "smiles": "C1C2CC3CC1CC(C2)C3", "mw": 136, "logp": 3.9, "tpsa": 0},
            {"name": "Cubane-rigid", "smiles": "C12C3C4C1C5C2C3C45", "mw": 104, "logp": 2.8, "tpsa": 0},
        ],
        "A": [  # Aromatic
            {"name": "Benzene", "smiles": "c1ccccc1", "mw": 78, "logp": 2.1, "tpsa": 0},
            {"name": "Pyridine", "smiles": "c1ccncc1", "mw": 79, "logp": 0.7, "tpsa": 13},
            {"name": "Imidazole", "smiles": "c1c[nH]cn1", "mw": 68, "logp": -0.2, "tpsa": 29},
        ],
        "M": [  # Allosteric/Covalent
            {"name": "Allosteric-scaffold", "smiles": "c1ccc2c(c1)nc(n2)N", "mw": 145, "logp": 1.2, "tpsa": 65},
            {"name": "Covalent-warhead", "smiles": "C=CC(=O)N", "mw": 71, "logp": -0.8, "tpsa": 43},
            {"name": "Beta-lactam", "smiles": "O=C1NC(=O)C1", "mw": 85, "logp": -1.1, "tpsa": 58},
        ],
    }
    
    @classmethod
    def search(cls, group: str, max_results: int = 10) -> List[Dict]:
        """
        Fake search - returns predefined molecules for each group
        In real version: query ZINC database with RDKit filters
        """
        molecules = cls.FAKE_MOLECULES.get(group, [])
        
        # Add some random variations
        results = []
        for mol in molecules[:max_results]:
            # Create slight variations
            for i in range(3):
                variant = mol.copy()
                variant["id"] = f"FAKE{group}{random.randint(1000, 9999)}"
                variant["mw"] = variant["mw"] * random.uniform(0.9, 1.1)
                variant["logp"] = variant["logp"] * random.uniform(0.8, 1.2)
                variant["score"] = random.uniform(0.7, 0.95)
                results.append(variant)
        
        # Sort by score
        results.sort(key=lambda x: x["score"], reverse=True)
        return results[:max_results]
    
    @classmethod
    def validate_molecule(cls, smiles: str, group: str) -> Dict:
        """
        Fake validation - always returns semi-random result
        In real version: use RDKit to calculate descriptors
        """
        # Simple SMILES length heuristic (totally fake!)
        fake_mw = len(smiles) * 12
        fake_logp = (smiles.count('c') - smiles.count('O') - smiles.count('N')) * 0.5
        fake_tpsa = (smiles.count('O') + smiles.count('N')) * 20
        
        violations = []
        
        # Fake Lipinski checks
        if fake_mw > 500:
            violations.append(f"MW {fake_mw:.0f} > 500")
        if fake_logp > 5:
            violations.append(f"LogP {fake_logp:.1f} > 5")
        
        # Fake group-specific checks
        criteria = cls._get_criteria(group)
        if "logp_min" in criteria and fake_logp < criteria["logp_min"]:
            violations.append(f"LogP {fake_logp:.1f} < {criteria['logp_min']}")
        if "logp_max" in criteria and fake_logp > criteria["logp_max"]:
            violations.append(f"LogP {fake_logp:.1f} > {criteria['logp_max']}")
        
        return {
            "pass": len(violations) == 0,
            "scores": {
                "MW": round(fake_mw, 1),
                "LogP": round(fake_logp, 2),
                "TPSA": round(fake_tpsa, 1),
                "HBD": smiles.count('O'),
                "HBA": smiles.count('O') + smiles.count('N'),
            },
            "violations": violations,
            "note": "⚠️ DEMO: Using simplified SMILES heuristics, not real RDKit calculations"
        }
    
    @classmethod
    def _get_criteria(cls, group: str) -> Dict:
        """Get filter criteria for each group"""
        criteria_map = {
            "H": {"logp_min": 3.0, "logp_max": 5.0},
            "P": {"logp_min": 0.0, "logp_max": 2.5},
            "C+": {"logp_min": 0.0, "logp_max": 3.0},
            "C-": {"logp_min": 0.5, "logp_max": 3.0},
            "V": {"logp_min": 0.5, "logp_max": 3.0},
            "A": {"logp_min": 2.0, "logp_max": 4.5},
            "M": {"logp_min": 1.0, "logp_max": 4.0},
        }
        return criteria_map.get(group, {})


class FilterCriteriaHelper:
    """Human-readable descriptions of filter criteria"""
    
    CRITERIA = {
        "H": {
            "name": "Increased Hydrophobicity",
            "description": "Mutations that increase hydrophobicity (e.g., R→I, K→L)",
            "target": "Hydrophobic molecules to fill hydrophobic pockets",
            "properties": {
                "LogP": "3.0–5.0 (highly lipophilic)",
                "TPSA": "10–70 Ų (low polarity)",
                "MW": "250–450 Da",
                "Aromatic rings": "≥1 (π-π stacking)",
                "HBD/HBA": "0–2 / 0–5 (minimal H-bonding)"
            }
        },
        "P": {
            "name": "Increased Polarity",
            "description": "Mutations that increase polarity (e.g., I→S, F→T)",
            "target": "Polar molecules to exploit new H-bond opportunities",
            "properties": {
                "LogP": "0.0–2.5 (hydrophilic)",
                "TPSA": "70–120 Ų (high polarity)",
                "HBD": "≥2 (H-bond donors)",
                "HBA": "4–10 (H-bond acceptors)"
            }
        },
        "C+": {
            "name": "Positive Charge Increase",
            "description": "Mutations adding positive charge (e.g., D→K, E→R)",
            "target": "Anionic molecules (COO⁻, SO₃⁻) to restore salt bridges",
            "properties": {
                "Formal charge": "−1 or −2",
                "pKa (acidic)": "≤6.5 (deprotonated at pH 7.4)",
                "LogP": "0.0–3.0"
            }
        },
        "C-": {
            "name": "Negative Charge Increase",
            "description": "Mutations adding negative charge (e.g., K→E, R→D)",
            "target": "Cationic molecules (NH₃⁺, guanidinium) for electrostatic interaction",
            "properties": {
                "Formal charge": "+1",
                "pKa (basic)": "≥7.5 (protonated at pH 7.4)",
                "LogP": "0.5–3.0"
            }
        },
        "V": {
            "name": "Steric/Volume Change",
            "description": "Large volume changes (e.g., G→W, W→A)",
            "target": "Size-matched molecules (small if cavity smaller, bulky if larger)",
            "properties": {
                "If d_vol > 0.6": "MW ≤300, flexible (RotB ≥4)",
                "If d_vol < −0.6": "MW ≥300, bulky scaffolds"
            }
        },
        "A": {
            "name": "Aromatic Loss",
            "description": "Loss of aromatic residues (e.g., F→A, Y→S)",
            "target": "Aromatic molecules to restore π-π stacking",
            "properties": {
                "Aromatic rings": "≥1",
                "LogP": "2.0–4.5",
                "Planarity": "Preferred for stacking"
            }
        },
        "M": {
            "name": "Motif/Active-Site Core",
            "description": "Mutations in critical functional regions",
            "target": "Allosteric inhibitors or covalent warheads (use with caution)",
            "properties": {
                "Strategy": "Target distant sites or use reactive groups",
                "MW": "300–500 Da",
                "Selectivity": "Critical - high risk of toxicity"
            }
        }
    }
    
    @classmethod
    def get_description(cls, group: str) -> str:
        """Get formatted description for a group"""
        if group not in cls.CRITERIA:
            return f"Unknown group: {group}"
        
        info = cls.CRITERIA[group]
        lines = [
            f"**{info['name']}**",
            f"",
            f"📋 {info['description']}",
            f"",
            f"🎯 **Design Strategy:** {info['target']}",
            f"",
            f"⚗️ **Target Properties:**"
        ]
        
        for key, val in info['properties'].items():
            lines.append(f"  • {key}: {val}")
        
        return "\n".join(lines)


# Quick test
if __name__ == "__main__":
    # Test classification
    features = {
        "d_hydro": -3.5,
        "d_vol": 0.2,
        "d_charge": 0.0,
        "rel_pos": 0.45,
        "center_dist": 0.1,
        "motif_pen": 0.0
    }
    
    group = MutationGroupClassifier.classify(features, "R", "I")
    print(f"R→I classified as: {group}")
    
    # Test fake molecule search
    molecules = FakeMoleculeDatabase.search(group, max_results=5)
    print(f"\nFound {len(molecules)} molecules for group {group}:")
    for mol in molecules:
        print(f"  {mol['id']}: {mol['name']} (MW={mol['mw']:.1f}, LogP={mol['logp']:.2f}, Score={mol['score']:.2f})")
    
    # Test validation
    test_smiles = "c1ccccc1"
    result = FakeMoleculeDatabase.validate_molecule(test_smiles, group)
    print(f"\nValidation of '{test_smiles}' for group {group}:")
    print(f"  Pass: {result['pass']}")
    print(f"  Scores: {result['scores']}")
    if result['violations']:
        print(f"  Violations: {', '.join(result['violations'])}")