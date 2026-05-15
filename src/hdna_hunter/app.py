#!/usr/bin/env python3
"""
app.py — Streamlit web UI for hdna-hunter.

Run:
    streamlit run src/hdna_hunter/app.py

This is a self-contained Streamlit application.  The heavy lifting is
performed by the hdna_hunter Python API (which calls the compiled C
extension internally).  No subprocess calls or external binaries are
required.
"""

from __future__ import annotations

import io
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st
import plotly.express as px

import hdna_hunter

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(page_title="H-DNA Hunter", page_icon="🧬", layout="wide")

st.markdown("""
<style>
    .main-title  { font-size: 2.4rem; font-weight: 700; color: #1a6b3c; }
    .subtitle    { font-size: 1.1rem; color: #555; margin-bottom: 1.5rem; }
    .section     { font-size: 1.15rem; font-weight: 600; color: #1a6b3c; margin-top: 1rem; }
    div[data-testid="stMetric"] {
        background: #f0f7f4; border-radius: 10px; padding: 12px;
    }
    .stButton > button {
        background-color: #1a6b3c; color: white;
        font-size: 1.1rem; font-weight: 600;
        padding: 0.6rem 2.5rem; border-radius: 8px; border: none;
    }
    .stButton > button:hover { background-color: #145c32; }
</style>
""", unsafe_allow_html=True)

st.markdown('<p class="main-title">🧬 H-DNA Hunter</p>', unsafe_allow_html=True)
st.markdown(
    '<p class="subtitle">H-DNA / Triplex Mirror Repeat Detector — '
    f'imperfect repeat support · v{hdna_hunter.__version__}</p>',
    unsafe_allow_html=True,
)

# ---------------------------------------------------------------------------
# Demo sequence (chr1 region with known H-DNA motifs)
# ---------------------------------------------------------------------------

DEMO_SEQ = """\
>chr1:43585222-43586222
GGAACTGCGTTCCTTTGGAGGAGGAGAGGCGCTCTGCGTTTTAGAGTTTCCAGTTTTTCT
GTTCTGTTTTTTCCCCATCTTTGTGGTTTTATCAACTTTTGGTCTTTGATGATGGTGATG
TACAGATGGGTTTTTGGTGTGGATGTCCTTTCTGTTTGTTAGTTTTCCATCTAACAGACA
GGACCCTCAGCTGCAAGTCTGTTGGAATACCCTGCTGTGTGAGGTGTCAGTGTGCCCCTG
CTGGGGGGTGCCTCCCAGTTAGGCTGCTCAGGGGTCAGGGGTCAGGGACCCACTTGAGGA
GGCAGTCTGCCCGTTCCCAGATCTCCAGCTGCGTGCTGGAAGAACCACTACTCTCTTCAA
AGCTGTCAGACAGGGACATTTAAGTCTGCAGAGGTTACTGCTGTCTTTTTGTTTGTCTGT
GCCCTGCCCCCAGAGGTGGAGCCTACAGAGGCAGGCAGGCCTCCTTGAGCTGTGGTGGGT
TCCACCCAGTTGGAGCTTCCCGGCTGCTTTGTTTACCTAAGCAAGCCTGGGCAATGGCGG
GCGCCCCTCCCCCAGCCTCCCTGCCGCCTTGCAGTTTGATCTCAGACTGCTGTGCTAGCA
ATCAGCGAGACTCCGTGGGCGTAGGACCCTCCGAGCCAGGTGCGGGATATCATCTCGTGG
TGCGCCGTTTAAGCCAGTCGGAAAAGCGCAGTATTCGGGTGGGAGTGACCCGATTTTCCA
GGTGCGTCCGTCACCCCTTTCTTTGACTCGGAAAGGGAACTCCCTGACCCCTTGCGCTTC
CCAAGTGAGGCAGTGCCTTGCCCTGCTTCGGCTCGTGCACGGTGCGCGCACCCACTGACC
TGCGCCCACTGTCTGGCACTCCCTAGTGAGATGAACCCGGTACCTCAGATGGAAATGCAG
AAATCACCCGTCTTCTGCGTCGCTCACGCTGGGAGCTGTAGATTCGGCCATCTTGGCTCC
TCCCCGACTTTGTTTTTCGATTGTTTTGTCCCTTGAGATT
"""

# ---------------------------------------------------------------------------
# Sidebar: parameters
# ---------------------------------------------------------------------------

with st.sidebar:
    st.markdown("## ⚙️ Parameters")

    minrep    = st.number_input("Min arm length (minrep)",  min_value=4,   max_value=100, value=10, step=1)
    maxrep    = st.number_input("Max arm length (maxrep)",  min_value=4,   max_value=200, value=50, step=1)
    maxspacer = st.number_input("Max spacer length",        min_value=0,   max_value=50,  value=7,  step=1)
    purity    = st.number_input(
        "Purity threshold", min_value=0.0, max_value=1.0, value=0.80,
        step=0.05, format="%.2f",
        help="Min fraction of GA or CT bases in each arm. Use 1.0 for strict mode.",
    )
    mismatch  = st.number_input(
        "Mismatch tolerance", min_value=0.0, max_value=1.0, value=0.10,
        step=0.05, format="%.2f",
        help="Max fraction of mirror-position mismatches allowed. Use 0.0 for exact mirror only.",
    )

    st.markdown("---")
    st.markdown("### 🔬 Post-filter")
    at_threshold = st.number_input(
        "Max AT content in either arm", min_value=0.0, max_value=1.0,
        value=0.80, step=0.05, format="%.2f",
        help="Drop hits where AT% in the left OR right arm exceeds this value.",
    )

    st.markdown("---")
    skip_overlap = st.checkbox("Skip overlap removal", value=False)
    st.markdown("**Strict** = purity 1.0 · mismatch 0.0  \n**Relaxed** = purity 0.8 · mismatch 0.2")

# ---------------------------------------------------------------------------
# Input section
# ---------------------------------------------------------------------------

st.markdown('<p class="section">📂 Input sequence</p>', unsafe_allow_html=True)
st.caption("Upload a FASTA file **or** paste below. Uploaded file takes precedence.")

col_up, col_paste = st.columns([1, 2])

with col_up:
    st.markdown("**Upload FASTA**")
    uploaded = st.file_uploader(
        "FASTA file", type=["fa", "fna", "fasta", "txt"],
        label_visibility="collapsed",
    )
    if uploaded:
        st.success(f"📄 {uploaded.name} ({uploaded.size:,} bytes)")

with col_paste:
    st.markdown("**Paste sequence**")
    pasted = st.text_area(
        "FASTA text", value=DEMO_SEQ, height=180,
        label_visibility="collapsed",
        help="Must start with a >header line.",
    )

if uploaded is not None:
    fasta_text = uploaded.getvalue().decode("utf-8")
    input_label = f"file: {uploaded.name}"
elif pasted and pasted.strip():
    fasta_text = pasted
    input_label = "pasted sequence"
else:
    fasta_text = None
    input_label = None

# ---------------------------------------------------------------------------
# Run button
# ---------------------------------------------------------------------------

st.markdown("")
run = st.button("🔬 Submit")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _at_pct(seq: str) -> float:
    s = seq.upper().replace(".", "")
    return (sum(1 for b in s if b in "AT") / len(s) * 100) if s else 0.0


def _run_scan(fasta_text: str) -> pd.DataFrame:
    """Write text to a temp file and call scan_fasta()."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
        fh.write(fasta_text)
        tmp_path = fh.name

    hits = hdna_hunter.scan_fasta(
        tmp_path,
        minrep=minrep,
        maxrep=maxrep,
        maxspacer=maxspacer,
        purity=purity,
        mismatch=mismatch,
        remove_overlaps=not skip_overlap,
    )

    Path(tmp_path).unlink(missing_ok=True)
    return pd.DataFrame(hits) if hits else pd.DataFrame()


# ---------------------------------------------------------------------------
# Run logic — results stored in session_state to survive widget interactions
# ---------------------------------------------------------------------------

if run:
    if not fasta_text:
        st.warning("Please upload a FASTA file or paste a sequence first.")
        st.stop()

    with st.spinner("Scanning for H-DNA motifs…"):
        try:
            df_raw = _run_scan(fasta_text)
        except Exception as exc:
            st.error(f"Scan failed: {exc}")
            st.stop()

    if df_raw.empty:
        st.session_state["results"] = None
        st.info("No hits found.")
        st.stop()

    # AT filter
    df_raw["left_at_pct"]  = df_raw["left_arm"].apply(_at_pct)
    df_raw["right_at_pct"] = df_raw["right_arm"].apply(_at_pct)
    n_before = len(df_raw)
    df_raw = df_raw[
        (df_raw["left_at_pct"]  <= at_threshold * 100) &
        (df_raw["right_at_pct"] <= at_threshold * 100)
    ].reset_index(drop=True)
    n_dropped = n_before - len(df_raw)

    st.session_state["results"]     = df_raw
    st.session_state["n_dropped"]   = n_dropped
    st.session_state["input_label"] = input_label


# ---------------------------------------------------------------------------
# Display results
# ---------------------------------------------------------------------------

if "results" in st.session_state and st.session_state["results"] is not None:
    df       = st.session_state["results"]
    n_hits   = len(df)
    n_dropped = st.session_state["n_dropped"]

    st.caption(f"Input source: **{st.session_state['input_label']}**")
    if n_dropped > 0:
        st.info(f"**{n_dropped:,}** hit(s) removed by AT filter. **{n_hits:,}** remaining.")

    # -- metrics row
    col_a, col_b, col_c, col_d = st.columns(4)
    col_a.metric("Total hits",        n_hits)
    col_b.metric("Unique seq IDs",    df["seq_id"].nunique() if "seq_id" in df.columns else 1)
    col_c.metric("Perfect mirrors",   int(df["is_perfect"].sum()))
    col_d.metric("Median arm (bp)",   int(df["arm_length"].median()))

    # -- results table
    st.markdown("---")
    st.markdown('<p class="section">🗂️ Results</p>', unsafe_allow_html=True)

    display_cols = [
        "seq_id", "start", "end",
        "arm_length", "spacer_length",
        "ga_pct", "ct_pct", "left_at_pct", "right_at_pct",
        "mirror_identity", "is_perfect",
        "left_arm", "spacer", "right_arm",
    ]
    display_cols = [c for c in display_cols if c in df.columns]
    float_cols   = ["ga_pct", "ct_pct", "mirror_identity", "left_at_pct", "right_at_pct"]

    display_df = df[display_cols].head(50).copy()
    for c in float_cols:
        if c in display_df.columns:
            display_df[c] = display_df[c].round(2)

    st.dataframe(display_df, use_container_width=True, height=420)
    st.caption(f"Showing first 50 of {n_hits:,} hits.")

    tsv_bytes = df.to_csv(sep="\t", index=False).encode()
    st.download_button(
        "⬇️ Download full TSV",
        data=tsv_bytes,
        file_name="hdna_results.tsv",
        mime="text/tab-separated-values",
    )

    # -- plots
    st.markdown("---")
    show_plots = st.checkbox("📈 Show distribution plots", value=False)

    if show_plots:
        col1, col2, col3 = st.columns(3)

        with col1:
            fig = px.histogram(
                df, x="arm_length", nbins=20, opacity=0.75,
                title="Arm length distribution",
                labels={"arm_length": "Arm length (bp)"},
                color_discrete_sequence=["#1a6b3c"],
            )
            fig.update_layout(margin=dict(t=40, b=20))
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            spacer_counts = (
                df.groupby("spacer_length").size().reset_index(name="count")
            )
            fig2 = px.bar(
                spacer_counts, x="spacer_length", y="count",
                title="Spacer length distribution",
                labels={"spacer_length": "Spacer (bp)", "count": "Count"},
                color_discrete_sequence=["#1a6b3c"],
            )
            fig2.update_layout(margin=dict(t=40, b=20))
            st.plotly_chart(fig2, use_container_width=True)

        with col3:
            fig3 = px.histogram(
                df, x="mirror_identity", nbins=20, opacity=0.75,
                title="Mirror identity distribution",
                labels={"mirror_identity": "Mirror identity (%)"},
                color_discrete_sequence=["#2196a6"],
            )
            fig3.update_layout(margin=dict(t=40, b=20))
            st.plotly_chart(fig3, use_container_width=True)

        col4, col5 = st.columns(2)

        with col4:
            fig4 = px.scatter(
                df, x="arm_length", y="mirror_identity",
                color="spacer_length",
                opacity=0.6,
                title="Arm length vs mirror identity",
                labels={
                    "arm_length": "Arm length (bp)",
                    "mirror_identity": "Mirror identity (%)",
                    "spacer_length": "Spacer",
                },
            )
            fig4.update_layout(margin=dict(t=40, b=20))
            st.plotly_chart(fig4, use_container_width=True)

        with col5:
            perfect_counts = df["is_perfect"].value_counts().reset_index()
            perfect_counts.columns = ["is_perfect", "count"]
            perfect_counts["label"] = perfect_counts["is_perfect"].map(
                {True: "Perfect mirror", False: "Imperfect mirror"}
            )
            fig5 = px.pie(
                perfect_counts, values="count", names="label",
                title="Perfect vs imperfect mirrors",
                color_discrete_sequence=["#1a6b3c", "#d1d5db"],
                hole=0.3,
            )
            fig5.update_layout(margin=dict(t=40, b=20))
            st.plotly_chart(fig5, use_container_width=True)
