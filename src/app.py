import streamlit as st
import pandas as pd
import altair as alt

# --- Import Modules ---
from sequence_fetcher import SequenceFetcher
from gc_content import calculate_gc_content
from kmer_frequency import calculate_kmer_frequency, get_most_common_kmers
from data_loader import load_gene_metadata, search_genes
from analysis_advanced import translate_dna, calculate_shannon_entropy, rolling_entropy, find_motif
from analysis_protein import calculate_protein_properties
from analysis_dna import calculate_gc_skew, calculate_codon_usage

# --- New Modules for Phase 2 ---
from languages import TRANSLATIONS
from disease_data import get_disease_variants

# --- Page Config ---
st.set_page_config(
    page_title="Human Genome Analysis",
    layout="wide",
    initial_sidebar_state="expanded",
    page_icon="🧬"
)

# --- Caching ---
@st.cache_data(show_spinner=False)
@st.cache_data(show_spinner=False)
def get_gene_metadata_v2():
    try:
        df = load_gene_metadata()
        if df is None:
            return "Function returned None (Check console for 'Error: File not found' or 'Error loading data')"
        return df
    except Exception as e:
        return str(e)

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

# --- Helper for Translation ---
def get_text(key, lang_code):
    return TRANSLATIONS.get(lang_code, TRANSLATIONS['EN']).get(key, key)

# --- Sidebar ---
st.sidebar.markdown("### 🌐 " + get_text('language_select', 'EN'))
lang_selection = st.sidebar.selectbox(
    "Select Language", 
    ["English (EN)", "Italiano (IT)", "Slovenský (SK)"], 
    label_visibility="collapsed"
)

# Extract language code
if "Italiano" in lang_selection:
    LANG = 'IT'
elif "Slovenský" in lang_selection:
    LANG = 'SK'
else:
    LANG = 'EN'

# Shortcut function for current language
def t(key):
    return get_text(key, LANG)

# --- Main Layout ---

# Title
col_title, col_logo = st.columns([4, 1])
with col_title:
    st.title(t('title'))
    st.markdown(t('subtitle'))

st.sidebar.header(t('sidebar_header'))

search_mode = st.sidebar.radio(t('search_method'), [t('search_by_name'), t('accession_id')], horizontal=True)

selected_id = ""

if search_mode == t('search_by_name'):
    df_genes = get_gene_metadata_v2()
    # Check if it's a DataFrame
    if isinstance(df_genes, pd.DataFrame):
        query = st.sidebar.text_input(t('search_by_name'), placeholder=t('gene_name_placeholder'))
        if query:
            results = search_genes(df_genes, query)
            if not results.empty:
                # Keep technical ID visible but maybe localized description if we had it
                options = results.apply(lambda x: f"{x['ID']} | {x['Description'][:40]}...", axis=1).tolist()
                selected_option = st.sidebar.selectbox(t('results'), options)
                selected_id = selected_option.split(" | ")[0]
            else:
                st.sidebar.warning(t('no_matches'))
    else:
        # It's an error message
        st.sidebar.error(f"{t('metadata_error')}\n\nError: {df_genes}")
        from data_loader import DEFAULT_DATA_PATH
        st.sidebar.caption(f"Checked path: `{DEFAULT_DATA_PATH}`")
else:
    selected_id = st.sidebar.text_input(t('accession_id'), value="NM_000014.6")

st.sidebar.markdown("---")
st.sidebar.header(t('settings_header'))
k_mer_len = st.sidebar.slider(t('kmer_len'), 1, 6, 3)
window_size = st.sidebar.slider(t('window_size'), 50, 500, 100, step=50)

fetch_btn = st.sidebar.button(t('analyze_btn'), type="primary")

# --- Main Content ---
if fetch_btn and selected_id:
    with st.spinner(f"{t('fetching')} {selected_id}..."):
        try:
            record = fetch_sequence_cached(selected_id)
        except Exception:
            record = None

    if record:
        # Summary Box
        st.success(f"{t('loaded')}: **{record.id}**")
        with st.expander(f"📝 {t('description')} & {t('raw_seq')}", expanded=False):
            st.markdown(f"**{t('description')}:** {record.description}")
            st.text_area("Sequence", str(record.seq)[:500] + "...", height=100)

        # Tabs
        tabs = st.tabs([
            t('tab_overview'), 
            t('tab_dna'), 
            t('tab_kmers'), 
            t('tab_protein'), 
            t('tab_advanced'),
            t('tab_disease') # New Tab
        ])

        # Tab 1: Overview
        with tabs[0]:
            col1, col2 = st.columns(2)
            with col1:
                st.metric(t('seq_len'), f"{len(record.seq):,}")
            
            gc = calculate_gc_content(record.seq)
            with col2:
                st.metric(t('gc_content'), f"{gc:.2f}%")
            
            # Simple Composition Chart
            bases = {'A': record.seq.count('A'), 'T': record.seq.count('T'), 
                     'G': record.seq.count('G'), 'C': record.seq.count('C')}
            df_bases = pd.DataFrame(list(bases.items()), columns=['Base', 'Count'])
            
            chart_base = alt.Chart(df_bases).mark_arc(innerRadius=50).encode(
                theta='Count',
                color=alt.Color('Base', scale=alt.Scale(domain=['A', 'T', 'G', 'C'], range=['#FF9AA2', '#B5EAD7', '#FFDAC1', '#C7CEEA'])),
                tooltip=['Base', 'Count']
            ).properties(title=t('base_comp'))
            st.altair_chart(chart_base, use_container_width=True)

        # Tab 2: DNA Structure
        with tabs[1]:
            st.subheader(t('gc_skew'))
            st.caption(t('gc_skew_desc'))
            
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
            st.subheader(t('codon_usage'))
            codon_data = get_codon_usage(record.seq)
            df_codons = pd.DataFrame(codon_data)
            
            with st.expander(t('codon_table')):
                st.dataframe(df_codons, use_container_width=True)

        # Tab 3: K-mers
        with tabs[2]:
            st.subheader(f"{k_mer_len}{t('kmer_freq')}")
            kmers = calculate_kmer_frequency(record.seq, k_mer_len)
            top_kmers = get_most_common_kmers(kmers, 15)
            
            if top_kmers:
                df_kmers = pd.DataFrame(top_kmers, columns=['K-mer', 'Count'])
                chart_kmers = alt.Chart(df_kmers).mark_bar().encode(
                    x=alt.X('K-mer', sort='-y'),
                    y='Count',
                    color='Count',
                    tooltip=['K-mer', 'Count']
                )
                st.altair_chart(chart_kmers, use_container_width=True)
            else:
                st.info(t('no_kmers'))

        # Tab 4: Protein Analysis
        with tabs[3]:
            st.subheader(t('trans_props'))
            protein_seq = translate_dna(record.seq)
            
            st.text_area(t('prot_seq'), str(protein_seq), height=150)
            
            if len(protein_seq) > 0:
                props = calculate_protein_properties(protein_seq)
                
                p_col1, p_col2, p_col3, p_col4 = st.columns(4)
                p_col1.metric(t('mol_weight'), f"{props.get('Molecular Weight',0):.0f} Da")
                p_col2.metric(t('isoelectric'), f"{props.get('Isoelectric Point',0):.2f}")
                p_col3.metric(t('instability'), f"{props.get('Instability Index',0):.2f}")
                p_col4.metric(t('aromaticity'), f"{props.get('Aromaticity',0):.2f}")
                
                if props.get('Instability Index', 0) > 40:
                    st.warning(t('unstable'))
                else:
                    st.success(t('stable'))

        # Tab 5: Advanced
        with tabs[4]:
            col_ent, col_mot = st.columns(2)
            
            with col_ent:
                st.subheader(t('entropy'))
                entropy = calculate_shannon_entropy(record.seq)
                st.metric(t('global_entropy'), f"{entropy:.4f} bits")
                
                if st.button(t('calc_rolling')):
                    df_ent = rolling_entropy(record.seq, window_size, 20)
                    chart_ent = alt.Chart(df_ent).mark_line(color='#FFA07A').encode(
                        x='Position',
                        y='Entropy'
                    )
                    st.altair_chart(chart_ent, use_container_width=True)

            with col_mot:
                st.subheader(t('motif_search'))
                motif_q = st.text_input(t('find_motif'), "TATA")
                if motif_q:
                    locs = find_motif(record.seq, motif_q)
                    if locs:
                        st.success(t('found_matches').format(len(locs)))
                        st.write(f"Positions: {locs[:10]}...")
                    else:
                        st.warning(t('not_found'))
                        
        # Tab 6: Disease Associations (New)
        with tabs[5]:
            st.subheader(t('disease_map'))
            st.info(t('disease_desc'))
            
            variants = get_disease_variants(record.id)
            
            if variants:
                df_variants = pd.DataFrame(variants)
                
                # Display table
                st.dataframe(
                    df_variants.rename(columns={
                        'pos': t('position'),
                        'ref': 'Ref',
                        'alt': 'Alt',
                        'disease': t('disease'),
                        'sig': t('significance'),
                        'desc': t('description')
                    }),
                    use_container_width=True
                )
                
                # Visualization of variants on the sequence timeline
                base_chart = alt.Chart(df_variants).mark_circle(size=100).encode(
                    x=alt.X('pos', title=t('position')),
                    y=alt.Y('disease', title=t('disease')),
                    color=alt.Color('sig', legend=alt.Legend(title=t('significance')), scale=alt.Scale(scheme='reds')),
                    tooltip=['pos', 'ref', 'alt', 'disease', 'desc']
                ).properties(height=200)
                
                st.altair_chart(base_chart, use_container_width=True)
                
            else:
                st.warning(f"No mock disease data available for gene {record.id}. Try searching for 'Insulin' (NM_000014).")

    else:
        st.error(f"{t('error_fetch')} {selected_id}")

elif not fetch_btn:
    st.info(t('landing_info'))
    
    st.markdown(t('features_title'))
    st.markdown(t('features_list'))
