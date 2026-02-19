import streamlit as st
import pandas as pd
import altair as alt
from sequence_fetcher import SequenceFetcher
from gc_content import calculate_gc_content
from kmer_frequency import calculate_kmer_frequency, get_most_common_kmers
from data_loader import load_gene_metadata, search_genes
from analysis_advanced import translate_dna, calculate_shannon_entropy, rolling_entropy, find_motif
from analysis_protein import calculate_protein_properties
from analysis_dna import calculate_gc_skew, calculate_codon_usage

# --- Page Config ---
st.set_page_config(
    page_title="Human Genome Analysis",
    layout="wide",
    initial_sidebar_state="expanded",
    page_icon="🧬"
)

# --- Caching ---
@st.cache_data
def get_gene_metadata():
    return load_gene_metadata()

@st.cache_data
def fetch_sequence_cached(accession_id):
    fetcher = SequenceFetcher()
    return fetcher.fetch_sequence(accession_id)

@st.cache_data
def get_gc_skew(seq, window, step):
    return calculate_gc_skew(seq, window, step)

@st.cache_data
def get_codon_usage(seq):
    return calculate_codon_usage(seq)

# --- Main Layout ---

# Title Area
col_title, col_logo = st.columns([4, 1])
with col_title:
    st.title("🧬 Human Genome Analysis")
    st.markdown("_Explore and analyze DNA sequences with advanced bioinformatics tools._")

# --- Sidebar ---
st.sidebar.header("🔍 Gene Selection")

search_mode = st.sidebar.radio("Search Method", ["Search by Name", "Accession ID"], horizontal=True)

selected_id = ""

if search_mode == "Search by Name":
    df_genes = get_gene_metadata()
    if df_genes is not None:
        query = st.sidebar.text_input("Gene Name", placeholder="e.g. Insulin, TP53")
        if query:
            results = search_genes(df_genes, query)
            if not results.empty:
                options = results.apply(lambda x: f"{x['ID']} | {x['Description'][:40]}...", axis=1).tolist()
                selected_option = st.sidebar.selectbox("Results", options)
                selected_id = selected_option.split(" | ")[0]
            else:
                st.sidebar.warning("No matches found.")
    else:
        st.sidebar.error("Metadata not available.")
else:
    selected_id = st.sidebar.text_input("Accession ID", value="NM_000014.6")

st.sidebar.markdown("---")
st.sidebar.header("⚙️ Analysis Settings")
k_mer_len = st.sidebar.slider("K-mer Length", 1, 6, 3)
window_size = st.sidebar.slider("Rolling Window (bp)", 50, 500, 100, step=50)

fetch_btn = st.sidebar.button("🚀 Analyze Sequence", type="primary")

# --- Main Content ---
if fetch_btn and selected_id:
    with st.spinner(f"Fetching sequence data for {selected_id}..."):
        record = fetch_sequence_cached(selected_id)

    if record:
        # Summary Box
        st.success(f"Loaded: **{record.id}**")
        with st.expander("📝 Sequence Description & Raw Data", expanded=True):
            st.markdown(f"**Description:** {record.description}")
            st.text_area("Raw Sequence (First 500bp)", str(record.seq)[:500] + "...", height=100)

        # Tabs
        tabs = st.tabs([
            "📊 Overview", 
            "🧬 DNA Structure", 
            "🦠 K-mers", 
            "🥩 Protein Analysis", 
            "🔎 Motifs & Entropy"
        ])

        # Tab 1: Overview
        with tabs[0]:
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Sequence Length", f"{len(record.seq):,} bp")
            
            gc = calculate_gc_content(record.seq)
            with col2:
                st.metric("GC Content", f"{gc:.2f}%")
            
            # Simple Composition Chart
            bases = {'A': record.seq.count('A'), 'T': record.seq.count('T'), 
                     'G': record.seq.count('G'), 'C': record.seq.count('C')}
            df_bases = pd.DataFrame(list(bases.items()), columns=['Base', 'Count'])
            
            chart_base = alt.Chart(df_bases).mark_arc(innerRadius=50).encode(
                theta='Count',
                color=alt.Color('Base', scale=alt.Scale(domain=['A', 'T', 'G', 'C'], range=['#FF9AA2', '#B5EAD7', '#FFDAC1', '#C7CEEA'])),
                tooltip=['Base', 'Count']
            ).properties(title="Base Composition")
            st.altair_chart(chart_base, use_container_width=True)

        # Tab 2: DNA Structure (GC Skew & Codon Usage)
        with tabs[1]:
            st.subheader("GC Skew Analysis")
            st.caption(f"Rolling GC Skew (Window={window_size}) to identify replication origins or biases.")
            
            skew_data = get_gc_skew(record.seq, window_size, step=20)
            df_skew = pd.DataFrame(skew_data)
            
            if not df_skew.empty:
                chart_skew = alt.Chart(df_skew).mark_area(
                    line={'color':'teal'},
                    color=alt.Gradient(
                        gradient='linear',
                        stops=[alt.GradientStop(color='white', offset=0),
                               alt.GradientStop(color='teal', offset=1)],
                        x1=1, x2=1, y1=1, y2=0
                    )
                ).encode(
                    x=alt.X('Position', title='Position (bp)'),
                    y=alt.Y('GC_Skew', title='GC Skew'),
                    tooltip=['Position', 'GC_Skew']
                ).properties(height=300)
                st.altair_chart(chart_skew, use_container_width=True)
            
            st.divider()
            st.subheader("Codon Usage Bias")
            codon_data = get_codon_usage(record.seq)
            df_codons = pd.DataFrame(codon_data)
            
            with st.expander("View Codon Usage Table"):
                st.dataframe(df_codons, use_container_width=True)

        # Tab 3: K-mers
        with tabs[2]:
            st.subheader(f"{k_mer_len}-mer Frequency Analysis")
            kmers = calculate_kmer_frequency(record.seq, k_mer_len)
            top_kmers = get_most_common_kmers(kmers, 15)
            
            if top_kmers:
                df_kmers = pd.DataFrame(top_kmers, columns=['K-mer', 'Count'])
                chart_kmers = alt.Chart(df_kmers).mark_bar().encode(
                    x=alt.X('K-mer', sort='-y'),
                    y='Count',
                    color=alt.Color('Count', scale=alt.Scale(scheme='tealblues')), # Changed scheme
                    tooltip=['K-mer', 'Count']
                )
                st.altair_chart(chart_kmers, use_container_width=True)
            else:
                st.info("No K-mers found.")

        # Tab 4: Protein Analysis
        with tabs[3]:
            st.subheader("Translation & Physicochemical Properties")
            protein_seq = translate_dna(record.seq)
            
            st.text_area("Protein Sequence", str(protein_seq), height=150)
            
            if len(protein_seq) > 0:
                props = calculate_protein_properties(protein_seq)
                
                p_col1, p_col2, p_col3, p_col4 = st.columns(4)
                p_col1.metric("Mol. Weight", f"{props.get('Molecular Weight',0):.0f} Da")
                p_col2.metric("Isoelectric Pt (pI)", f"{props.get('Isoelectric Point',0):.2f}")
                p_col3.metric("Instability Idx", f"{props.get('Instability Index',0):.2f}")
                p_col4.metric("Aromaticity", f"{props.get('Aromaticity',0):.2f}")
                
                if props.get('Instability Index', 0) > 40:
                    st.warning("⚠️ Protein is classified as **unstable** (Index > 40).")
                else:
                    st.success("✅ Protein is classified as **stable**.")

        # Tab 5: Advanced (Entropy & Motifs)
        with tabs[4]:
            col_ent, col_mot = st.columns(2)
            
            with col_ent:
                st.subheader("Entropy Analysis")
                entropy = calculate_shannon_entropy(record.seq)
                st.metric("Global Shannon Entropy", f"{entropy:.4f} bits")
                
                if st.button("Calculate Rolling Entropy"):
                    df_ent = rolling_entropy(record.seq, window_size, 20)
                    chart_ent = alt.Chart(df_ent).mark_line(color='#FFA07A').encode(
                        x='Position',
                        y='Entropy'
                    )
                    st.altair_chart(chart_ent, use_container_width=True)

            with col_mot:
                st.subheader("Motif Search")
                motif_q = st.text_input("Find Motif (DNA)", "TATA")
                if motif_q:
                    locs = find_motif(record.seq, motif_q)
                    if locs:
                        st.success(f"Found {len(locs)} matches.")
                        st.write(f"First 10 positions: {locs[:10]}")
                    else:
                        st.warning("Not found.")

    else:
        st.error(f"Could not fetch data for ID: {selected_id}")

elif not fetch_btn:
    # Landing Page Info
    st.info("👈 Please select a gene or enter an Accession ID in the sidebar to start.")
    
    st.markdown("""
    ### Features
    - **Protein Analysis**: Molecular weight, pI, and stability prediction.
    - **DNA Structure**: GC Skew visualization for detecting replication origins.
    - **Interactive Plots**: Zoomable charts for K-mers and entropy.
    """)
