from fastapi import FastAPI
from pydantic import BaseModel, Field
from typing import Optional
from .pipeline import mutational_scan

app = FastAPI(title="TB Resistance Predictor (MVP)")

class PredictIn(BaseModel):
    sequence: str = Field(..., description="Protein sequence, 1-letter amino acids")
    region_start: int = Field(..., description="Start position (1-indexed, AA position)")
    region_end: int = Field(..., description="End position (1-indexed, AA position)")
    gene_hint: Optional[str] = Field(None, description="e.g., rpoB, katG, pncA ...")
    nucleotide_sequence: Optional[str] = Field(None, description="Optional: original DNA/RNA sequence for NT-level mutations")
    orf_start: int = Field(0, description="ORF start position (0-indexed) in nucleotide sequence")
    enable_drug_filtering: bool = Field(True, description="Enable 6-group drug filtering")
    resistance_cutoff: float = Field(0.7, description="Resistance probability cutoff")

@app.post("/predict")
def predict(payload: PredictIn):
    res = mutational_scan(
        sequence=payload.sequence,
        region_start=payload.region_start,
        region_end=payload.region_end,
        gene_hint=payload.gene_hint or "",
        enable_drug_filtering=payload.enable_drug_filtering,
        resistance_cutoff=payload.resistance_cutoff,
        demo_mode=True,
        nucleotide_sequence=payload.nucleotide_sequence,
        orf_start=payload.orf_start
    )
    return res

@app.get("/health")
def health():
    return {"ok": True}
