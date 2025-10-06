import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

from tb_resistance_mvp.pipeline import mutational_scan
try:
    from tb_resistance_mvp.drug_filtering_demo import FakeMoleculeDatabase, FilterCriteriaHelper
    DEMO_AVAILABLE = True
except ImportError:
    DEMO_AVAILABLE = False

# Import translation functions (inline to avoid new file dependency)
def translate_dna_sequence(sequence: str, frame: int = 0) -> str:
    """Translate DNA sequence to protein"""
    CODON_TABLE = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    
    sequence = sequence.upper().replace(" ", "").replace("\n", "")
    sequence = sequence.replace("U", "T")
    
    valid_bases = set("ATCG")
    if not all(base in valid_bases for base in sequence):
        raise ValueError("Invalid DNA sequence. Only A, T, C, G allowed.")
    
    sequence = sequence[frame:]
    
    protein = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    
    return ''.join(protein)


def detect_orf(sequence: str, min_length: int = 100):
    """Detect Open Reading Frames"""
    sequence = sequence.upper().replace("U", "T").replace(" ", "").replace("\n", "")
    orfs = []
    
    for frame in range(3):
        seq_frame = sequence[frame:]
        
        for i in range(0, len(seq_frame) - 2, 3):
            codon = seq_frame[i:i+3]
            
            if codon == 'ATG':
                start = i + frame
                
                for j in range(i, len(seq_frame) - 2, 3):
                    stop_codon = seq_frame[j:j+3]
                    if stop_codon in ('TAA', 'TAG', 'TGA'):
                        end = j + frame + 3
                        
                        if (end - start) >= min_length:
                            orf_seq = sequence[start:end]
                            try:
                                protein = translate_dna_sequence(orf_seq, frame=0)
                                orfs.append({
                                    'start': start + 1,
                                    'end': end,
                                    'frame': frame,
                                    'length_nt': end - start,
                                    'length_aa': len(protein),
                                    'protein': protein
                                })
                            except:
                                pass
                        break
    
    return orfs


st.set_page_config(page_title="TB Resistance Predictor — Demo v2.0", layout="wide", initial_sidebar_state="expanded")

# Custom CSS
st.markdown("""
<style>
    .demo-badge {
        background-color: #ff6b6b;
        color: white;
        padding: 4px 12px;
        border-radius: 12px;
        font-size: 12px;
        font-weight: bold;
    }
    .molecule-card {
        border: 1px solid #444;
        border-radius: 8px;
        padding: 12px;
        margin: 8px 0;
        background-color: #1e1e1e;
        color: #fafafa;
    }
    .molecule-card code {
        background-color: #2d2d2d;
        color: #61dafb;
        padding: 2px 6px;
        border-radius: 4px;
    }
    .molecule-card small {
        color: #b0b0b0;
    }
</style>
""", unsafe_allow_html=True)

# Header
col_h1, col_h2 = st.columns([4, 1])
with col_h1:
    st.title("🧬 TB Resistance Predictor — Demo v2.0")
with col_h2:
    if DEMO_AVAILABLE:
        st.markdown('<span class="demo-badge">DEMO MODE</span>', unsafe_allow_html=True)

st.caption("**NEW:** 6-group mutation classification & drug candidate suggestions (using simulated data for demo)")

# Sidebar
with st.sidebar:
    st.header("ℹ️ About Demo Mode")
    st.info("""
    **Demo Features:**
    - ✅ Full 6-group classification
    - ✅ Filter criteria display
    - ✅ Example molecule suggestions
    - ✅ Nucleotide → Protein translation
    - ⚠️ Using simulated SMILES & properties
    
    **For Production:**
    Replace with real RDKit + ZINC database
    """)
    
    st.header("📚 Mutation Groups")
    with st.expander("Group Descriptions"):
        st.markdown("""
        - **H**: Hydrophobic increase
        - **P**: Polar increase  
        - **C+**: Positive charge gain
        - **C-**: Negative charge gain
        - **V**: Volume change
        - **A**: Aromatic loss
        - **M**: Motif/Core region
        """)

# Session state
if "res" not in st.session_state:
    st.session_state.res = None
if "cutoff" not in st.session_state:
    st.session_state.cutoff = 0.70
if "show_adv" not in st.session_state:
    st.session_state.show_adv = False
if "enable_filtering" not in st.session_state:
    st.session_state.enable_filtering = True

# --- INPUTS ---
st.markdown("### 📝 Input Configuration")

# NEW: Input type selector
input_type = st.radio(
    "Input sequence type:",
    ["Amino Acid (Protein)", "Nucleotide (DNA/RNA)"],
    horizontal=True,
    help="Choose whether to input protein or nucleotide sequence"
)

protein_seq = ""

if input_type == "Amino Acid (Protein)":
    # Original AA input
    seq = st.text_area(
        "Protein sequence (1-letter AA)", 
        height=120, 
        value="MASTKQLLAVGHVPRNTLDEYQFGLSWVKMTPNIAERCLVDSGHTVQFPARLKMGYDETIVNQALSRPWFKGCYMDTNSLEIRAQPLHGMSFYTKVDACIGLPNRHTQVEWMDASFKRLYTPGQNVLSEACRKIFHDGTMVLPQANSRWDYEKGTMLFPISHVQRDANCLGTEYWMKPRSFAVLQNDGSITRKMLEHVPV",
        help="Enter your protein sequence using single-letter amino acid codes"
    )
    protein_seq = seq.strip().upper()
    
else:  # Nucleotide input
    nt_seq = st.text_area(
        "Nucleotide sequence (DNA/RNA)", 
        height=120, 
        value="ATGGCGTCAACAAAGCAGCTGCTGGCGGTCCATGTGCCGCGTAACACTGGACGAATATCAGTTTGGCCTGAGCTGGGTGAAGATGACGCCGAACATCGCCGAGCGCTGCCTGGTCGATTCGGGCCACACGGTGCAGTTCCCGGCGCGGCTCAAGATGGGCTACGATGAAACCATCGTGAACCAGGCGCTGAGCCGTCCGTGGTTCAAGGGCTGCTACATGGACACCAACAGCCTGGAGATCCGGGCGCAGCCGCTGCACGGGATGTCGTTCTACACGAAGGTGGACGCCTGCATCGGCCTGCCGAACCGGCACACGCAGGTGGAGTGGATGGACGCCTCATTCAAGCGGCTGTACACCCCGGGCCAGAACGTGCTGTCGGAGGCCTGCCGGAAGATCTTCCACGACGGCACCATGGTGCTGCCGCAGGCGAACAGCCGCTGGGACTATGAGAAGGGCACGATGCTGTTCCCGATCAGCCACGTCCAGCGCGACGCGAACTGCCTGGGCACCGAATACTGGATGAAGCCGCGCAGCTTCGCCGTGCTGCAGAACGACGGCAGCATCACCCGCAAGATGCTGGAGCACGTGCCGGTG",
        help="Enter DNA or RNA sequence (will be translated to protein)"
    )
    
    # Translation options
    col_frame, col_auto = st.columns(2)
    
    with col_frame:
        frame = st.selectbox(
            "Reading frame",
            [0, 1, 2],
            help="Select reading frame for translation (0 = start from position 1)"
        )
    
    with col_auto:
        auto_detect = st.checkbox(
            "Auto-detect ORF",
            value=True,
            help="Automatically detect the longest Open Reading Frame"
        )
    
    # Translate
    try:
        if auto_detect:
            orfs = detect_orf(nt_seq, min_length=100)
            if orfs:
                longest_orf = max(orfs, key=lambda x: x['length_aa'])
                protein_seq = longest_orf['protein']
                st.success(f"✅ Detected ORF: {longest_orf['length_aa']} AA (Frame {longest_orf['frame']}, Pos {longest_orf['start']}-{longest_orf['end']})")
                
                # Show all ORFs if multiple found
                if len(orfs) > 1:
                    with st.expander(f"🔍 Found {len(orfs)} ORFs - Click to view all"):
                        orf_data = []
                        for orf in orfs:
                            orf_data.append({
                                "Frame": orf['frame'],
                                "Start": orf['start'],
                                "End": orf['end'],
                                "Length (nt)": orf['length_nt'],
                                "Length (AA)": orf['length_aa']
                            })
                        st.dataframe(pd.DataFrame(orf_data), use_container_width=True)
            else:
                st.error("❌ No ORF detected (min length: 100 nt). Try manual frame selection.")
                protein_seq = ""
        else:
            protein_seq = translate_dna_sequence(nt_seq, frame=frame)
            st.info(f"ℹ️ Translated: {len(protein_seq)} amino acids (Frame {frame})")
        
        # Show translated sequence
        if protein_seq:
            with st.expander("🔬 View Translated Protein Sequence"):
                st.code(protein_seq, language="text")
                st.caption(f"Length: {len(protein_seq)} amino acids")
                
    except ValueError as e:
        st.error(f"❌ Translation error: {e}")
        protein_seq = ""
    except Exception as e:
        st.error(f"❌ Unexpected error: {e}")
        protein_seq = ""

# Region selection
col1, col2, col3 = st.columns(3)
with col1:
    if protein_seq:
        max_pos = len(protein_seq)
        default_start = min(50, max_pos)
        default_end = min(100, max_pos)
    else:
        max_pos = 1000
        default_start = 50
        default_end = 100
    
    start = st.number_input("Region start (1-indexed)", min_value=1, max_value=max_pos, value=default_start, step=1)
with col2:
    end = st.number_input("Region end (1-indexed)", min_value=1, max_value=max_pos, value=default_end, step=1)
with col3:
    gene = st.text_input("Gene hint (optional)", value="rpoB", placeholder="e.g., rpoB, katG, pncA")

# Settings
with st.expander("⚙️ Analysis Settings", expanded=False):
    col_a, col_b = st.columns(2)
    with col_a:
        st.session_state.cutoff = st.slider(
            "Resistance probability cutoff",
            min_value=0.0, max_value=1.0, step=0.05, value=float(st.session_state.cutoff),
            help="Mutations with P(resist) ≥ cutoff are classified as High risk"
        )
    with col_b:
        st.session_state.enable_filtering = st.checkbox(
            "Enable drug filtering (6-group)",
            value=bool(st.session_state.enable_filtering),
            help="Classify mutations and suggest drug candidates"
        )
    
    st.session_state.show_adv = st.toggle(
        "Show advanced physicochemical columns",
        value=bool(st.session_state.show_adv)
    )

scan = st.button("🔬 Scan Mutations", type="primary", use_container_width=True)


def build_df(res, cutoff):
    """Build DataFrame from results"""
    rows = []
    for r in res["matrix"]:
        base = {
            "pos": r["pos"], 
            "wt": r["wt"], 
            "mut": r["mut"], 
            "mechanism": r["mechanism"],
            "group": r.get("group", "N/A")
        }
        base.update(r["features"])
        for d, p in r["p_resist"].items():
            base[f"P({d})"] = round(p, 4)
        rows.append(base)
    
    df = pd.DataFrame(rows)
    p_cols = [c for c in df.columns if re.match(r'^P\(.+\)$', c)]
    
    if p_cols:
        df["Probability (max)"] = df[p_cols].max(axis=1).round(4)
        df["Risk"] = np.where(df["Probability (max)"] >= cutoff, "High", "Low")
        
        def drugs_above_cutoff(row):
            lst = [d for d in p_cols if row[d] >= cutoff]
            return ", ".join(d.replace("P(", "").replace(")", "") for d in lst) if lst else ""
        df["Drugs ≥ cutoff"] = df.apply(drugs_above_cutoff, axis=1)
    else:
        df["Probability (max)"] = 0.0
        df["Risk"] = "Low"
        df["Drugs ≥ cutoff"] = ""
    
    return df, p_cols


def render_molecule_card(mol, group):
    """Render molecule card with dark background"""
    st.markdown(f"""
    <div class="molecule-card">
        <strong>{mol.get('name', 'Unknown')}</strong> 
        <span style="float:right; color:#888;">Score: {mol.get('score', 0):.2f}</span>
        <br>
        <code>{mol.get('smiles', 'N/A')}</code>
        <br>
        <small>
        MW: {mol.get('mw', 0):.1f} | 
        LogP: {mol.get('logp', 0):.2f} | 
        TPSA: {mol.get('tpsa', 0):.1f}
        {f" | Charge: {mol.get('charge', 0):+d}" if 'charge' in mol else ""}
        </small>
    </div>
    """, unsafe_allow_html=True)


def render_results(res, cutoff, show_adv, enable_filtering):
    """Render complete results"""
    st.markdown("---")
    st.markdown("### 📊 Mutation Scan Results")
    
    df, p_cols = build_df(res, cutoff)
    
    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total mutations", len(df))
    with col2:
        high_risk = (df["Risk"] == "High").sum()
        pct = high_risk/len(df)*100 if len(df) > 0 else 0
        st.metric("High risk", high_risk, delta=f"{pct:.1f}%")
    with col3:
        if enable_filtering and "drug_filtering" in res:
            num_groups = len(res["drug_filtering"]["filter_summary"])
            st.metric("Groups identified", num_groups)
        else:
            st.metric("Groups", "Disabled")
    with col4:
        avg_prob = df["Probability (max)"].mean()
        st.metric("Avg P(resist)", f"{avg_prob:.3f}")
    
    # Mutation table
    st.markdown("#### 📋 Mutation Table")
    
    advanced_cols = ["d_hydro", "d_vol", "d_charge", "rel_pos", "center_dist", "motif_pen"]
    base_cols = ["pos", "wt", "mut", "group", "mechanism", "Probability (max)", "Risk", "Drugs ≥ cutoff"] + p_cols
    visible_cols = base_cols + (advanced_cols if show_adv else [])
    visible_cols = [c for c in visible_cols if c in df.columns]
    
    st.dataframe(df[visible_cols], use_container_width=True, height=420)
    
    # Heatmaps
    st.markdown("### 📈 Resistance Escape Maps")
    
    num_drugs = len(res["drugs"])
    cols_per_row = 3
    
    for row_idx in range((num_drugs + cols_per_row - 1) // cols_per_row):
        cols = st.columns(cols_per_row)
        for col_idx in range(cols_per_row):
            drug_idx = row_idx * cols_per_row + col_idx
            if drug_idx >= num_drugs:
                break
            
            d = res["drugs"][drug_idx]
            hm = res["heatmap"][d]
            xs = list(hm.keys())
            ys = list(hm.values())
            
            with cols[col_idx]:
                fig, ax = plt.subplots(figsize=(5, 3))
                ax.plot(xs, ys, linewidth=2, color='steelblue')
                ax.axhline(y=cutoff, color='red', linestyle='--', alpha=0.5, label=f'Cutoff ({cutoff})')
                ax.fill_between(xs, ys, cutoff, where=[y >= cutoff for y in ys], alpha=0.3, color='red')
                ax.set_xlabel("Position", fontsize=9)
                ax.set_ylabel(f"Max P({d})", fontsize=9)
                ax.set_title(f"{d}", fontsize=11, fontweight='bold')
                ax.legend(fontsize=7)
                ax.grid(alpha=0.3)
                st.pyplot(fig)
                plt.close()
    
    # Drug filtering
    if enable_filtering and "drug_filtering" in res and res["drug_filtering"]["high_risk_count"] > 0:
        st.markdown("---")
        st.markdown("### 🧪 Drug Candidate Filtering Strategy")
        
        df_summary = res["drug_filtering"]
        
        if df_summary.get("demo_mode"):
            st.warning("⚠️ **DEMO MODE**: Using simulated molecules. Replace with real ZINC database for production.")
        
        st.success(f"🎯 **{df_summary['high_risk_count']} high-risk mutations** identified (≥{cutoff:.2f} probability)")
        
        if df_summary["filter_summary"]:
            tabs = st.tabs([f"Group {g}" for g in df_summary["filter_summary"].keys()])
            
            for idx, (group, info) in enumerate(df_summary["filter_summary"].items()):
                with tabs[idx]:
                    st.markdown(f"### {info['name']}")
                    st.caption(f"**{info['mutation_count']} mutations** in this group")
                    
                    st.markdown("#### 📋 Filter Criteria")
                    st.markdown(info["description"])
                    
                    st.markdown("#### 🧬 Example Mutations")
                    mut_cols = st.columns(5)
                    for i, mut in enumerate(info["example_mutations"][:5]):
                        with mut_cols[i % 5]:
                            st.code(mut)
                    
                    st.markdown("---")
                    st.markdown("#### 💊 Example Drug Candidates (from simulated DB)")
                    
                    if "example_molecules" in df_summary:
                        molecules = df_summary["example_molecules"].get(group, [])
                        
                        if molecules:
                            for mol in molecules[:5]:
                                render_molecule_card(mol, group)
                        else:
                            st.info("No example molecules available for this group")
    
    # Download
    st.markdown("---")
    st.markdown("### 📥 Export Results")
    
    col_down1, col_down2 = st.columns(2)
    
    with col_down1:
        df_full, _ = build_df(res, cutoff)
        csv = df_full.to_csv(index=False).encode("utf-8")
        st.download_button(
            "📄 Download Full Mutation Table (CSV)",
            csv,
            file_name="mutation_scan_results.csv",
            mime="text/csv",
            use_container_width=True
        )
    
    with col_down2:
        if enable_filtering and "drug_filtering" in res:
            candidate_data = []
            for group, mols in df_summary.get("example_molecules", {}).items():
                for mol in mols:
                    candidate_data.append({
                        "group": group,
                        "molecule_id": mol.get("id", ""),
                        "name": mol.get("name", ""),
                        "smiles": mol.get("smiles", ""),
                        "mw": mol.get("mw", 0),
                        "logp": mol.get("logp", 0),
                        "score": mol.get("score", 0)
                    })
            
            if candidate_data:
                candidates_csv = pd.DataFrame(candidate_data).to_csv(index=False).encode("utf-8")
                st.download_button(
                    "💊 Download Drug Candidates (CSV)",
                    candidates_csv,
                    file_name="drug_candidates.csv",
                    mime="text/csv",
                    use_container_width=True
                )


# Run scan
if scan:
    if not protein_seq:
        st.error("❌ Please provide a valid sequence (protein or nucleotide)")
    elif start > end:
        st.error("❌ Region start must be ≤ region end")
    elif end > len(protein_seq):
        st.error(f"❌ Region end ({end}) exceeds sequence length ({len(protein_seq)})")
    else:
        try:
            with st.spinner("🔍 Scanning mutations..."):
                st.session_state.res = mutational_scan(
                    protein_seq,
                    int(start),
                    int(end),
                    gene,
                    enable_drug_filtering=st.session_state.enable_filtering,
                    resistance_cutoff=st.session_state.cutoff,
                    demo_mode=True
                )
            st.success("✅ Scan completed successfully!")
            st.balloons()
        except Exception as e:
            st.error(f"❌ Error during scanning: {e}")
            st.exception(e)

# Display results
if st.session_state.res:
    render_results(
        st.session_state.res,
        st.session_state.cutoff,
        st.session_state.show_adv,
        st.session_state.enable_filtering
    )
else:
    st.info("👆 Configure your analysis settings above and click **Scan Mutations** to begin.")
    
    with st.expander("📖 Quick Tutorial", expanded=False):
        st.markdown("""
        ### How to use this tool:
        
        **New: Nucleotide Input Support! 🧬**
        - Choose "Nucleotide (DNA/RNA)" to input DNA sequence
        - Auto-detect ORF or manually select reading frame
        - System will translate to protein automatically
        
        **Standard Workflow:**
        1. **Select input type**: Amino Acid or Nucleotide
        2. **Enter sequence**: Paste your sequence
        3. **Set region**: Choose which positions to scan
        4. **Gene hint**: Specify gene name (optional)
        5. **Adjust settings**: Cutoff & drug filtering
        6. **Click Scan**: Run the analysis
        7. **Review results**: Mutations, groups, drug candidates
        8. **Download**: Export CSV for further analysis
        
        ### Understanding Groups:
        - **H (Hydrophobic)**: Need lipophilic molecules
        - **P (Polar)**: Need H-bond donors/acceptors
        - **C+ (Charge+)**: Need anionic molecules (COO⁻)
        - **C- (Charge-)**: Need cationic molecules (NH₃⁺)
        - **V (Volume)**: Size-matched scaffolds
        - **A (Aromatic)**: π-system molecules
        - **M (Motif)**: Allosteric or covalent inhibitors
        """)