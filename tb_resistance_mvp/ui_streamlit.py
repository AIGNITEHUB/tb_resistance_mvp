import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from tb_resistance_mvp.pipeline import mutational_scan

st.set_page_config(page_title="TB Resistance Predictor — MVP", layout="wide")
st.title("🧬 TB Resistance Predictor — MVP")

seq = st.text_area("Protein sequence (1-letter AA)", height=160, value="MASTKQLLAVGHVPRNTLDEYQFGLSWVKMTPNIAERCLVDSGHTVQFPARLKMGYDETIVNQALSRPWFKGCYMDTNSLEIRAQPLHGMSFYTKVDACIGLPNRHTQVEWMDASFKRLYTPGQNVLSEACRKIFHDGTMVLPQANSRWDYEKGTMLFPISHVQRDANCLGTEYWMKPRSFAVLQNDGSITRKMLEHVPV")
col1, col2, col3 = st.columns(3)
with col1:
    start = st.number_input("Region start (1-indexed)", min_value=1, value=200, step=1)
with col2:
    end = st.number_input("Region end (1-indexed)", min_value=1, value=300, step=1)
with col3:
    gene = st.text_input("Gene hint (optional, e.g., rpoB, pncA)", value="rpoB")

if st.button("Scan"):
    try:
        res = mutational_scan(seq, int(start), int(end), gene)
        # Build a DF for table view
        rows = []
        for r in res["matrix"]:
            base = {"pos": r["pos"], "wt": r["wt"], "mut": r["mut"], "mechanism": r["mechanism"]}
            base.update(r["features"])
            for d,p in r["p_resist"].items():
                base[f"P({d})"] = p
            rows.append(base)
        df = pd.DataFrame(rows)
        st.success(f"Scanned positions {res['region']['start']}–{res['region']['end']} (length={res['length']}) | Drugs: {', '.join(res['drugs'])}")
        st.dataframe(df, use_container_width=True, height=400)

        # Heatmap per-drug (matplotlib; one plot per drug, no seaborn, no styles)
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

        # --- NEW: Design suggestions ---
        st.subheader("🔬 Design suggestions (pilot)")

        # Regions table
        if "design_regions" in res and res["design_regions"]:
            dr_df = pd.DataFrame(res["design_regions"])
            st.markdown("**Vùng ưu tiên (escape thấp)**")
            st.dataframe(dr_df, use_container_width=True, height=180)
        else:
            st.info("Không tìm thấy vùng ưu tiên rõ rệt trong khoảng đã chọn.")

        # Scaffold & property hints
        if "design" in res:
            d = res["design"]
            c1, c2 = st.columns([2,1])
            with c1:
                st.markdown("**Seed scaffolds (SMILES)**")
                st.code("\n".join(f"- {s['name']}: {s['smiles']}  # {s['note']}" for s in d["seed_scaffolds"]))
            with c2:
                st.markdown("**Thuộc tính gợi ý (ước lượng)**")
                st.json(d["props"], expanded=True)


        # Offer CSV downloads
        csv = df.to_csv(index=False).encode("utf-8")
        st.download_button("Download mutation table (CSV)", csv, file_name="mut_table.csv", mime="text/csv")

    except Exception as e:
        st.error(f"Error: {e}")
