import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

from tb_resistance_mvp.pipeline import mutational_scan

st.set_page_config(page_title="TB Resistance Predictor — MVP ", layout="wide")
st.title("🧬 TB Resistance Predictor — MVP ")

# --- Session state init ---
if "res" not in st.session_state:
    st.session_state.res = None
if "cutoff" not in st.session_state:
    st.session_state.cutoff = 0.70
if "show_adv" not in st.session_state:
    st.session_state.show_adv = False

# --- Inputs ---
seq = st.text_area("Protein sequence (1-letter AA)", height=160, value="MASTKQLLAVGHVPRNTLDEYQFGLSWVKMTPNIAERCLVDSGHTVQFPARLKMGYDETIVNQALSRPWFKGCYMDTNSLEIRAQPLHGMSFYTKVDACIGLPNRHTQVEWMDASFKRLYTPGQNVLSEACRKIFHDGTMVLPQANSRWDYEKGTMLFPISHVQRDANCLGTEYWMKPRSFAVLQNDGSITRKMLEHVPV")
col1, col2, col3 = st.columns(3)
with col1:
    start = st.number_input("Region start (1-indexed)", min_value=1, value=150, step=1)
with col2:
    end = st.number_input("Region end (1-indexed)", min_value=1, value=200, step=1)
with col3:
    gene = st.text_input("Gene hint (optional, e.g., rpoB, pncA)", value="rpoB")

# --- Controls available BEFORE scanning ---
with st.expander("Display & classification settings", expanded=True):
    st.session_state.cutoff = st.number_input(
        "Cutoff (xác suất ≥ cutoff → High risk)",
        min_value=0.0, max_value=1.0, step=0.05, value=float(st.session_state.cutoff),
        help="Áp dụng cho cột Probability (max) và Risk."
    )
    st.session_state.show_adv = st.toggle(
        "Hiện các cột hoá–lý (d_hydro → motif_pen)",
        value=bool(st.session_state.show_adv)
    )

scan = st.button("Scan")

def build_df(res, cutoff):
    rows = []
    for r in res["matrix"]:
        base = {"pos": r["pos"], "wt": r["wt"], "mut": r["mut"], "mechanism": r["mechanism"]}
        base.update(r["features"])
        for d,p in r["p_resist"].items():
            base[f"P({d})"] = round(p,4)
        rows.append(base)
    df = pd.DataFrame(rows)
    p_cols = [c for c in df.columns if re.match(r'^P\(.+\)$', c)]
    if p_cols:
        df["Probability (max)"] = df[p_cols].max(axis=1).round(4)
        df["Risk"] = np.where(df["Probability (max)"] >= cutoff, "High", "Low")
        def drugs_above_cutoff(row):
            lst = [d for d in p_cols if row[d] >= cutoff]
            return ", ".join(d.replace("P(","").replace(")","") for d in lst) if lst else ""
        df["Drugs ≥ cutoff"] = df.apply(drugs_above_cutoff, axis=1)
    else:
        df["Probability (max)"] = 0.0
        df["Risk"] = "Low"
        df["Drugs ≥ cutoff"] = ""
    return df, p_cols

def render_results(res, cutoff, show_adv):
    st.markdown("### Kết quả quét đột biến")
    df, p_cols = build_df(res, cutoff)
    advanced_cols = ["d_hydro","d_vol","d_charge","rel_pos","center_dist","motif_pen"]
    base_cols = ["pos","wt","mut","mechanism","Probability (max)","Risk","Drugs ≥ cutoff"] + p_cols
    visible_cols = base_cols + (advanced_cols if show_adv else [])
    visible_cols = [c for c in visible_cols if c in df.columns]
    st.dataframe(df[visible_cols], use_container_width=True, height=420)

    # Heatmaps
    for d in res["drugs"]:
        hm = res["heatmap"][d]
        xs = list(hm.keys())
        ys = list(hm.values())
        fig, ax = plt.subplots()
        ax.plot(xs, ys)
        ax.set_xlabel("Position")
        ax.set_ylabel(f"Max P({d}) at pos")
        ax.set_title(f"Escape map — {d}")
        st.pyplot(fig)

    # Design suggestions
    st.subheader("🔬 Design suggestions (pilot)")
    if res.get("design_regions"):
        dr_df = pd.DataFrame(res["design_regions"])
        st.markdown("**Vùng ưu tiên (escape thấp)**")
        st.dataframe(dr_df, use_container_width=True, height=180)
    d = res.get("design", {})
    if d:
        c1, c2 = st.columns([2,1])
        with c1:
            st.markdown("**Seed scaffolds (SMILES)**")
            if "seed_scaffolds" in d:
                lines = [f"- {s['name']}: {s['smiles']}  # {s['note']}" for s in d["seed_scaffolds"]]
                st.code("\n".join(lines))
        with c2:
            st.markdown("**Thuộc tính gợi ý**")
            st.json(d.get("props", {}), expanded=True)

    # Download
    df_full, _ = build_df(res, cutoff)
    csv = df_full.to_csv(index=False).encode("utf-8")
    st.download_button("Download mutation table (CSV)", csv, file_name="mut_table.csv", mime="text/csv")

# --- Run scan if button pressed ---
if scan:
    try:
        st.session_state.res = mutational_scan(seq, int(start), int(end), gene)
    except Exception as e:
        st.error(f"Error: {e}")

# --- Persisted results show even when toggling or changing cutoff ---
if st.session_state.res:
    render_results(st.session_state.res, st.session_state.cutoff, st.session_state.show_adv)
else:
    st.info("Nhập chuỗi & vùng rồi bấm **Scan**. Bạn có thể chỉnh cutoff và bật/tắt cột nâng cao trước khi quét.")
