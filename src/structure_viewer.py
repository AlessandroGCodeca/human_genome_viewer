from stmol import showmol
import py3Dmol

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
