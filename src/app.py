import streamlit as st
import pandas as pd
import altair as alt
from sequence_fetcher import SequenceFetcher
from gc_content import calculate_gc_content
from kmer_frequency import calculate_kmer_frequency, get_most_common_kmers
from data_loader import load_gene_metadata, search_genes
from analysis_advanced import translate_dna, calculate_shannon_entropy, rolling_entropy, find_motif

# Set page config
st.set_page_config(page_title="Human Genome Analysis", layout="wide")

# --- Caching ---
@st.cache_data
def get_gene_metadata():
    return load_gene_metadata()

@st.cache_data
def fetch_sequence_cached(accession_id):
    fetcher = SequenceFetcher()
    return fetcher.fetch_sequence(accession_id)

# --- Title ---
st.title("Human Genome Analysis Dashboard")
st.markdown("Analyze molecular sequences from the human genome.")

# --- Sidebar ---
st.sidebar.header("Configuration")

# Search Mode Selection
search_mode = st.sidebar.radio("Input Method", ["Search by Gene Name", "Enter Accession ID"])

selected_id = ""

if search_mode == "Search by Gene Name":
    # Load metadata (cached)
    df_genes = get_gene_metadata()
    
    if df_genes is not None:
        query = st.sidebar.text_input("Search Gene (e.g. 'Insulin', 'BRCA1')", help="Case-insensitive search in descriptions")
        
        if query:
            results = search_genes(df_genes, query)
            if not results.empty:
                # Create detailed options for the selectbox
                options = results.apply(lambda x: f"{x['ID']} | {x['Description'][:50]}...", axis=1).tolist()
                selected_option = st.sidebar.selectbox("Select Result", options)
                # Extract ID from selection
                selected_id = selected_option.split(" | ")[0]
            else:
                st.sidebar.warning("No matches found.")
        else:
            st.sidebar.info("Enter a search term above.")
    else:
        st.sidebar.error("Could not load gene metadata.")
        
else: # Direct ID Input
    selected_id = st.sidebar.text_input("Accession ID", value="NM_000014.6", help="Enter a genbank accession ID (e.g., NM_...)")


fetch_btn = st.sidebar.button("Fetch and Analyze")

# --- Main Area ---
if fetch_btn and selected_id:
    with st.spinner(f"Fetching {selected_id}..."):
        # Use cached fetcher
        record = fetch_sequence_cached(selected_id)
        
    if record:
        st.success(f"Successfully fetched {selected_id}")
        st.markdown(f"**Description:** {record.description}")
        
        # Tabs for different analyses
        tab1, tab2, tab3, tab4 = st.tabs(["Overview", "K-mer Analysis", "Translation", "Advanced Analysis"])
        
        with tab1:
            st.header("Overview")
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Length (bp)", len(record.seq))
                
            # GC Content
            gc_content = calculate_gc_content(record.seq)
            with col2:
                st.metric("GC Content", f"{gc_content:.2f}%")
                st.progress(gc_content / 100.0)
                
            # Sequence Viewer
            st.subheader("Sequence Preview (First 1000 bp)")
            st.text_area("DNA Sequence", str(record.seq[:1000]), height=200)

        with tab2:
            st.header("K-mer Frequency Analysis")
            k_col1, k_col2 = st.columns([1, 3])
            
            with k_col1:
                k = st.slider("Select K-mer length", min_value=1, max_value=6, value=3)
                top_n = st.slider("Show top N K-mers", min_value=5, max_value=50, value=10)

            kmers = calculate_kmer_frequency(record.seq, k)
            top_kmers = get_most_common_kmers(kmers, top_n)
            
            with k_col2:
                if top_kmers:
                    df_kmers = pd.DataFrame(top_kmers, columns=['K-mer', 'Count'])
                    
                    # Altair Chart
                    chart = alt.Chart(df_kmers).mark_bar().encode(
                        x=alt.X('K-mer', sort='-y'),
                        y='Count',
                        tooltip=['K-mer', 'Count']
                    ).properties(title=f"Top {top_n} {k}-mers")
                    
                    st.altair_chart(chart, use_container_width=True)
                    
                    with st.expander("Show raw K-mer counts"):
                        st.write(df_kmers)
                else:
                    st.warning("No K-mers found.")

        with tab3:
            st.header("Sequence Translation (DNA -> Protein)")
            protein_seq = translate_dna(record.seq)
            
            st.metric("Protein Length (aa)", len(protein_seq))
            st.text_area("Protein Sequence", protein_seq, height=300)
            st.caption("*Note: Translation uses standard genetic code. '*' indicates stop codon.*")

        with tab4:
            st.header("Advanced Analysis")
            
            # Shannon Entropy
            st.subheader("Sequence Complexity (Shannon Entropy)")
            entropy = calculate_shannon_entropy(record.seq) # Corrected function name
            st.metric("Global Entropy (bits)", f"{entropy:.4f}")
            st.caption("Higher entropy indicates more randomness (complex sequence). Max for DNA is 2.0 bits.")
            
            st.subheader("Rolling Entropy")
            window_size = st.slider("Window Size", 50, 500, 100)
            step_size = st.slider("Step Size", 10, 100, 20)
            
            if st.button("Calculate Rolling Entropy"):
                df_entropy = rolling_entropy(record.seq, window_size, step_size)
                if not df_entropy.empty:
                    chart_ent = alt.Chart(df_entropy).mark_line().encode(
                        x='Position',
                        y=alt.Y('Entropy', scale=alt.Scale(domain=[0, 2])),
                        tooltip=['Position', 'Entropy']
                    ).properties(title=f"Rolling Entropy (Window={window_size})")
                    st.altair_chart(chart_ent, use_container_width=True)
                else:
                    st.warning("Sequence too short for this window size.")
            
            st.subheader("Motif Search")
            motif = st.text_input("Enter Motif to Search (e.g., GAATTC)", "GAATTC")
            if motif:
                indices = find_motif(record.seq, motif)
                if indices:
                    st.success(f"Found {len(indices)} occurrences of '{motif}'")
                    st.write(f"Positions (0-based): {indices[:50]} {'...' if len(indices) > 50 else ''}")
                else:
                    st.warning(f"Motif '{motif}' not found.")

    else:
        st.error(f"Failed to fetch sequence for ID: {selected_id}. Please check the ID and try again.")
elif fetch_btn and not selected_id:
    st.warning("Please select a sequence ID first.")
else:
    st.info("Use the sidebar to search for a gene or enter an ID to start analysis.")
    
    st.markdown("---")
    st.markdown("### About")
    st.markdown("This dashboard fetches sequences from NCBI and performs basic analysis.")
    st.markdown("**New:** Search by gene name, view translations, and analyze sequence entropy!")
