# TB Resistance Predictor — MVP

Minimal, offline-friendly MVP to scan a user-selected region of a 500-aa protein
and **simulate** resistance probabilities per drug for every single-amino-acid mutation.
This is a scaffold you can extend with real biological models (ESM/EVE/FoldX).

> ⚠️ This MVP uses simple, deterministic heuristics (physicochemical deltas) — **not** real resistance models.

## Quickstart

### 1) Create a venv & install deps
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### 2) Run the API
```bash
uvicorn tb_resistance_mvp.api:app --reload --port 8000
```
- Open: http://localhost:8000/docs

### 3) Run the Streamlit UI (optional)
```bash
streamlit run tb_resistance_mvp/ui_streamlit.py
```
- Paste a sequence, choose region (e.g., 200–300), click **Scan** to see a heatmap + table.

## Project layout

```
tb_resistance_mvp/
  __init__.py
  constants.py
  features.py       # physicochemical features per mutation
  model.py          # simple deterministic scorer (per-drug)
  pipeline.py       # end-to-end mutational scan
  api.py            # FastAPI: /predict
  ui_streamlit.py   # Streamlit app for quick demo
requirements.txt
README.md
```

## Replace the heuristic core with real models

- In `features.py`, add ESM/EVE/ΔΔG, pocket distances, conservation, etc.
- In `model.py`, swap the `HeuristicModel` with your trained XGBoost/MLP heads.
- Keep `pipeline.py` & `api.py` as-is — they only depend on a `.predict(features)` interface.

## Safety / Ethics
This demo is for **in-silico** exploration only. It does **not** provide wet-lab guidance.
