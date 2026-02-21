import py3Dmol
import requests
from typing import Optional, Tuple

def render_pdb(pdb_id, width=800, height=600):
    """
    Render a 3D protein structure from a PDB ID using py3Dmol and stmol.
    """
    view = py3Dmol.view(query=f'pdb:{pdb_id}', width=width, height=height)
    
    # Style: Cartoon representation colored by secondary structure
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    
    # Add surface (optional, maybe make interactive later)
    # view.addSurface(py3Dmol.VDW, {'opacity': 0.7, 'color': 'white'})
    
    view.zoomTo()
    view.spin(True) # Auto-spin for cool effect
    
    return view

def render_uploaded_pdb(pdb_string, width=800, height=600):
    """
    Render a PDB structure from a string content.
    """
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_string, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()
    view.spin(True)
    return view

def fetch_pdb_by_gene(gene_symbol: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Cross-references a gene symbol (e.g., TP53) to find an associated PDB ID
    using MyGene.info, then fetches the raw PDB text from RCSB.
    Returns (pdb_id, pdb_raw_text) or (None, None) if not found.
    """
    if not gene_symbol:
        return None, None

    try:
        # 1. Look up PDB ID from Gene Symbol or Accession
        mg_url = f"https://mygene.info/v3/query?q={gene_symbol}&species=human&fields=pdb"
        mg_res = requests.get(mg_url, timeout=5)
        mg_res.raise_for_status()
        hits = mg_res.json().get('hits', [])
        
        if not hits:
            return None, None
            
        # Find the first hit that actually contains PDB entries
        pdb_list = None
        for hit in hits:
            if 'pdb' in hit:
                pdb_list = hit['pdb']
                break
                
        if not pdb_list:
            return None, None
            
        # Extract the first PDB ID (often an array, sometimes a string if there's only one)
        pdb_id = pdb_list[0] if isinstance(pdb_list, list) else pdb_list
        pdb_id = str(pdb_id).lower()

        # 2. Fetch the actual PDB file content from RCSB
        rcsb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        rcsb_res = requests.get(rcsb_url, timeout=10)
        rcsb_res.raise_for_status()
        
        return pdb_id.upper(), rcsb_res.text

    except requests.RequestException as e:
        print(f"Error fetching PDB for {gene_symbol}: {e}")
        return None, None
