import requests
import pandas as pd
import streamlit as st
import altair as alt

@st.cache_data(show_spinner=False, ttl=86400)
def fetch_ensembl_transcripts(gene_symbol: str):
    """
    Query the Ensembl REST API for a specific human gene symbol
    to extract all associated transcripts (isoforms) and their exon structures.
    """
    server = "https://rest.ensembl.org"
    # Step 1: Lookup the gene symbol to get the Ensembl Gene ID
    ext = f"/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1"
    
    headers = {"Content-Type": "application/json"}
    try:
        r = requests.get(server + ext, headers=headers)
        if not r.ok:
            return None
            
        data = r.json()
        
        # Parse the transcripts
        transcripts = data.get('Transcript', [])
        
        exons_data = []
        for t in transcripts:
            t_id = t.get('id')
            t_name = t.get('display_name', t_id)
            biotype = t.get('biotype')
            
            # We mostly care about protein-coding isoforms for this viewer
            if biotype != 'protein_coding':
                continue
                
            for exon in t.get('Exon', []):
                exons_data.append({
                    "Transcript": t_name,
                    "Start": exon.get('start'),
                    "End": exon.get('end'),
                    "Strand": exon.get('strand')
                })
                
        if not exons_data:
            return None
            
        df = pd.DataFrame(exons_data)
        
        return df
        
    except Exception as e:
        print(f"Error fetching Isoforms: {e}")
        return None

def plot_isoforms(df: pd.DataFrame):
    """
    Renders an Altair Gantt chart showing transcript isoforms and their exons.
    """
    if df is None or df.empty:
        return None
        
    # Standardize start/end for plotting (Altair needs standard min/max orientation)
    df['PlotStart'] = df[['Start', 'End']].min(axis=1)
    df['PlotEnd'] = df[['Start', 'End']].max(axis=1)
    
    # Calculate global domain for X axis to ensure alignment
    min_x = df['PlotStart'].min()
    max_x = df['PlotEnd'].max()
    
    # Create the chart
    chart = alt.Chart(df).mark_bar(cornerRadius=3, height=15).encode(
        x=alt.X('PlotStart:Q', scale=alt.Scale(domain=[min_x, max_x]), title='Genomic Position (bp)'),
        x2='PlotEnd:Q',
        y=alt.Y('Transcript:N', sort='-x', title='', axis=alt.Axis(labels=True, ticks=False)),
        color=alt.Color('Transcript:N', legend=None, scale=alt.Scale(scheme='category20b')),
        tooltip=[
            alt.Tooltip('Transcript:N', title='Isoform'),
            alt.Tooltip('Start:Q', title='Exon Start', format=','),
            alt.Tooltip('End:Q', title='Exon End', format=',')
        ]
    ).properties(
        height=max(200, len(df['Transcript'].unique()) * 30), # Dynamic height
        title="Protein-Coding Isoforms (Exon Maps)"
    ).interactive(bind_y=False) # Allow zooming on X but keep Y fixed
    
    return chart
