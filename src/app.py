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
# --- Phase 2 & 3 Modules ---
from languages import TRANSLATIONS
from disease_data import get_disease_variants
from structure_viewer import render_pdb
from simulator import MutationSimulator
from alignment import perform_pairwise_alignment
from ai_chat import AIGeneAssistant # Phase 4

# --- Page Config ---
st.set_page_config(
    page_title="Human Genome Analysis",
    layout="wide",
    initial_sidebar_state="expanded",
    page_icon="🧬"
)

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
            t('tab_3d'),      
            t('tab_mutation'), 
            t('tab_align'),    
            t('tab_disease'),
            t('tab_ai'), # Phase 4
            t('tab_advanced')
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
            
            seq1_input = st.text_area(t('seq1'), str(record.seq)[:200], height=100)
            seq2_input = st.text_area(t('seq2'), str(record.seq)[:200], height=100)
            
            if st.button(t('align_btn')):
                res = perform_pairwise_alignment(
                    seq1_input.replace("\n",""), 
                    seq2_input.replace("\n","")
                )
                st.metric(t('alignment_score'), f"{res['score']:.1f}")
                
                st.text(f"Best Match:\n{res['alignment_str']}")

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
                
                # Population Genetics (Phase 4)
                st.divider()
                st.subheader(t('allele_freq'))
                
                # Select a variant to show frequencies for
                variant_options = [f"{v['disease']} ({v['ref']}>{v['alt']} @ {v['pos']})" for v in variants]
                selected_var_idx = st.selectbox(t('variant'), range(len(variants)), format_func=lambda x: variant_options[x], key='pop_gen_select')
                
                selected_var = variants[selected_var_idx]
                if 'frequencies' in selected_var:
                    freqs = selected_var['frequencies']
                    # Localize population names
                    pop_names = {
                        'AFR': t('pop_afr'), 'AMR': t('pop_amr'), 
                        'EAS': t('pop_eas'), 'EUR': t('pop_eur'), 
                        'SAS': t('pop_sas')
                    }
                    
                    df_freq = pd.DataFrame([
                        {'Population': pop_names.get(k, k), 'Frequency': v} 
                        for k, v in freqs.items()
                    ])
                    
                    chart_freq = alt.Chart(df_freq).mark_bar().encode(
                        x=alt.X('Population', sort='-y'),
                        y=alt.Y('Frequency'),
                        color=alt.Color('Population', legend=None),
                        tooltip=['Population', 'Frequency']
                    ).properties(height=300)
                    
                    st.altair_chart(chart_freq, use_container_width=True)
                else:
                    st.info("No population frequency data for this variant.")

        # Tab 7: AI Assistant (Phase 4)
        with tabs[6]:
            st.subheader(t('tab_ai'))
            st.write(t('ai_intro'))
            
            ai = AIGeneAssistant()
            
            # Initialize chat history
            if "messages" not in st.session_state:
                st.session_state.messages = []
                # Welcome message
                st.session_state.messages.append({"role": "assistant", "content": t('ai_welcome')})

            # Display chat messages
            for message in st.session_state.messages:
                with st.chat_message(message["role"]):
                    st.markdown(message["content"])

            # Chat input
            if prompt := st.chat_input(t('ai_placeholder')):
                # Add user message
                st.session_state.messages.append({"role": "user", "content": prompt})
                with st.chat_message("user"):
                    st.markdown(prompt)

                # Generate AI response
                response = ai.get_response(prompt, record.id)
                
                # Add assistant message
                st.session_state.messages.append({"role": "assistant", "content": response})
                with st.chat_message("assistant"):
                    st.markdown(response)

        # Tab 8: Advanced
        with tabs[7]: # Shifted index
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

    else:
        st.error(f"{t('error_fetch')} {selected_id}")

elif not fetch_btn:
    st.info(t('landing_info'))
    
    st.markdown(t('features_title'))
    st.markdown(t('features_list'))
