from fastapi import FastAPI
from pydantic import BaseModel, Field
from typing import Optional
from .pipeline import mutational_scan

app = FastAPI(title="TB Resistance Predictor (MVP)")

class PredictIn(BaseModel):
    sequence: str = Field(..., description="Protein sequence, 1-letter amino acids")
    region_start: int
    region_end: int
    gene_hint: Optional[str] = Field(None, description="e.g., rpoB, katG, pncA ...")

@app.post("/predict")
def predict(payload: PredictIn):
    res = mutational_scan(
        sequence=payload.sequence,
        region_start=payload.region_start,
        region_end=payload.region_end,
        gene_hint=payload.gene_hint or ""
    )
    return res

@app.get("/health")
def health():
    return {"ok": True}
