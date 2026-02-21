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

# --- Phase 2 & 3 Modules ---
from languages import TRANSLATIONS
from disease_data import get_disease_variants
from structure_viewer import render_pdb
from simulator import MutationSimulator
from alignment import perform_pairwise_alignment
from ai_assistant import get_ai_response

# --- Page Config ---
st.set_page_config(
    page_title="Human Genome Analysis",
    layout="wide",
    initial_sidebar_state="expanded",
    page_icon="🧬"
)

import os
def load_css():
    css_path = os.path.join(os.path.dirname(__file__), "assets", "style.css")
    if os.path.exists(css_path):
        with open(css_path) as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

# Inject custom CSS
load_css()

# --- Caching ---
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

if 'active_id' not in st.session_state:
    st.session_state.active_id = None

if fetch_btn and selected_id:
    st.session_state.active_id = selected_id

# --- Main Content ---
if st.session_state.active_id:
    active_id = st.session_state.active_id
    with st.spinner(f"{t('fetching')} {active_id}..."):
        try:
            record = fetch_sequence_cached(active_id)
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
            t('tab_3d'),      # New
            t('tab_mutation'), # New
            t('tab_align'),    # New
            t('tab_disease'),
            t('tab_advanced'),
            t('tab_ai')        # AI Assistant
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
            
            # Pass sequence as string to avoid Streamlit UnhashableParamError for Bio.Seq objects
            skew_data = get_gc_skew(str(record.seq), window_size, step=20)
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
            codon_data = get_codon_usage(str(record.seq))
            df_codons = pd.DataFrame(codon_data)
            
            with st.expander(t('codon_table')):
                st.dataframe(df_codons, use_container_width=True)

        # Tab 3: 3D Structure (New)
        with tabs[2]:
            st.subheader(t('tab_3d'))
            
            col_pdb, col_viewer = st.columns([1, 4])
            with col_pdb:
                pdb_id = st.text_input(t('pdb_input'), "1TRZ")
                render = st.button(t('render_3d'))
                
            with col_viewer:
                if render and pdb_id:
                    # from stmol import showmol
                    # view = render_pdb(pdb_id)
                    # showmol(view, height=500, width=800)
                    try:
                        from stmol import showmol
                        view = render_pdb(pdb_id)
                        showmol(view, height=500, width=800)
                    except ImportError:
                        st.error("Please install 'stmol' and 'py3Dmol' to view 3D structures.")
                    except Exception as e:
                         st.error(f"Error rendering PDB: {e}")

        # Tab 4: Mutation Lab (New)
        with tabs[3]:
            st.subheader(t('tab_mutation'))
            st.write(t('mutation_intro'))
            
            sim = MutationSimulator(record.seq)
            
            # Editable text area initialized with original sequence
            user_seq = st.text_area("Sequence Editor", str(record.seq), height=150)
            
            if st.button(t('simulate_btn')):
                sim.mutated_dna = user_seq.replace("\n", "").replace(" ", "").upper()
                res = sim.compare_protein_properties()
                
                c1, c2 = st.columns(2)
                with c1:
                    st.markdown(f"**{t('original')}**")
                    st.code(res['Protein_Original'])
                with c2:
                    st.markdown(f"**{t('mutated')}**")
                    st.code(res['Protein_Mutated'])
                
                if res['Is_Silent']:
                    st.info(t('silent_mutation'))
                elif res['Stop_Codon_Introduced']:
                    st.warning(t('stop_codon'))
                    
                st.subheader(t('protein_diff'))
                if res['Properties_Comparison']:
                    st.table(pd.DataFrame(res['Properties_Comparison']).T)
                else:
                    st.write("No significant property changes detected.")

        # Tab 5: Alignment (New)
        with tabs[4]:
            st.subheader(t('tab_align'))
            st.write(t('align_intro'))
            
            st.markdown(f"**{t('msa_input_desc')}**")
            
            # Default example sequences for the tree
            default_fasta = f">Selected_{record.id}\n{str(record.seq)}\n>Variant_1\n{str(record.seq)[:400] + 'TGAC' + str(record.seq)[404:]}\n>Variant_2\n{str(record.seq)[10:300] + 'AAATTT' + str(record.seq)[306:450]}\n>Distant_Variant\n{'A'*50 + str(record.seq)[50:200] + 'C'*50}"
            
            fasta_input = st.text_area("FASTA Input", default_fasta, height=300, label_visibility="collapsed")
            
            if st.button(t('generate_tree_btn')):
                # Parse simple FASTA
                fasta_dict = {}
                current_name = ""
                current_seq = []
                for line in fasta_input.split('\n'):
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith('>'):
                        if current_name:
                            fasta_dict[current_name] = "".join(current_seq)
                        current_name = line[1:]
                        current_seq = []
                    else:
                        current_seq.append(line.upper().replace(" ", ""))
                if current_name:
                    fasta_dict[current_name] = "".join(current_seq)
                
                if len(fasta_dict) < 3:
                    st.warning("Please provide at least 3 sequences for a meaningful phylogenetic tree.")
                else:
                    with st.spinner("Calculating distance matrix and building tree..."):
                        from alignment import build_phylogenetic_tree
                        fig, matrix_data = build_phylogenetic_tree(fasta_dict)
                        
                        st.subheader(t('phylogenetic_tree'))
                        st.pyplot(fig)
                        
                        with st.expander(t('distance_matrix')):
                            st.dataframe(pd.DataFrame(matrix_data), use_container_width=True)

        # Tab 6: Disease Associations
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
                        'desc': t('description'),
                        'freq_desc': t('freq_desc')
                    }).drop(columns=['allele_freq'], errors='ignore'),
                    use_container_width=True
                )
                
                # Visualization of variants on the sequence timeline
                base_chart = alt.Chart(df_variants).mark_circle(size=100).encode(
                    x=alt.X('pos', title=t('position')),
                    y=alt.Y('disease', title=t('disease')),
                    color=alt.Color('sig', legend=alt.Legend(title=t('significance')), scale=alt.Scale(scheme='reds')),
                    tooltip=['pos', 'ref', 'alt', 'disease', 'desc', 'freq_desc']
                ).properties(height=200)
                
                st.altair_chart(base_chart, use_container_width=True)
                
                # Population Genetics (Allele Frequencies)
                st.markdown(f"#### {t('allele_freq')}")
                
                # Flatten the allele frequency data for Altair
                freq_data = []
                for _, row in df_variants.iterrows():
                    if 'allele_freq' in row and isinstance(row['allele_freq'], dict):
                        for pop, freq in row['allele_freq'].items():
                            pop_label = t(f'pop_{pop.lower()}') if f'pop_{pop.lower()}' in TRANSLATIONS[LANG] else pop
                            freq_data.append({
                                'Variant': f"{row['disease']} ({row['pos']}{row['ref']}>{row['alt']})",
                                'Population': pop_label,
                                'Frequency': freq
                            })
                
                if freq_data:
                    df_freqs = pd.DataFrame(freq_data)
                    freq_chart = alt.Chart(df_freqs).mark_bar().encode(
                        x=alt.X('Population:N', title='Population', axis=alt.Axis(labelAngle=-45)),
                        y=alt.Y('Frequency:Q', title='Frequency'),
                        color='Population:N',
                        column=alt.Column('Variant:N', title='Variants'),
                        tooltip=['Variant', 'Population', 'Frequency']
                    ).properties(width=120, height=250)
                    
                    st.altair_chart(freq_chart)
                
            else:
                st.warning(f"No mock disease data available for gene {record.id}. Try searching for 'Insulin' (NM_000014).")

        # Tab 7: Advanced
        with tabs[6]: # Shifted index
            col_ent, col_mot = st.columns(2)
            # ... (Existing Advanced content)
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

        # Tab 8: AI Assistant
        with tabs[7]:
            st.subheader(t('tab_ai'))
            
            gemini_api_key = st.text_input("Enter your Google Gemini API Key:", type="password", key="gemini_key_input")
            st.caption("You can get a free API key from [Google AI Studio](https://aistudio.google.com/).")
            
            st.divider()
            
            # Use session state to store chat history for the current selected ID
            if "messages" not in st.session_state:
                st.session_state.messages = []
            
            # Reset messages if a newly searched gene ID is found, or we can just keep them.
            # To be simple, we just append to the state.
            
            # Display chat messages from history on app rerun
            for message in st.session_state.messages:
                with st.chat_message(message["role"]):
                    st.markdown(message["content"])
            
            # If empty, add a greeting message from the assistant
            if not st.session_state.messages:
                greeting = t('ai_greeting').format(record.id)
                st.session_state.messages.append({"role": "assistant", "content": greeting})
                with st.chat_message("assistant"):
                    st.markdown(greeting)

            # Accept user input
            if prompt := st.chat_input(t('ai_placeholder')):
                # Display user message in chat message container
                with st.chat_message("user"):
                    st.markdown(prompt)
                # Add user message to chat history
                st.session_state.messages.append({"role": "user", "content": prompt})
                
                # Display assistant response in chat message container
                with st.chat_message("assistant"):
                    with st.spinner("Thinking..."):
                        response = get_ai_response(active_id, prompt, api_key=gemini_api_key)
                    st.markdown(response)
                # Add assistant response to chat history
                st.session_state.messages.append({"role": "assistant", "content": response})

    else:
        st.error(f"{t('error_fetch')} {active_id}")

elif not st.session_state.active_id:
    st.info(t('landing_info'))
    
    st.markdown(t('features_title'))
    st.markdown(t('features_list'))
