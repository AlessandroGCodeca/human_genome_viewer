from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq

def calculate_protein_properties(protein_sequence):
    """
    Calculate physicochemical properties of a protein sequence.
    Returns a dictionary of properties.
    """
    # Clean sequence (remove stop codons if present as '*')
    clean_seq = str(protein_sequence).replace('*', '')
    
    if not clean_seq:
        return {
            'Molecular Weight': 0,
            'Isoelectric Point': 0,
            'Instability Index': 0,
            'Aromaticity': 0,
            'Gravy': 0
        }
    
    analysed_seq = ProteinAnalysis(clean_seq)
    
    try:
        props = {
            'Molecular Weight': analysed_seq.molecular_weight(),
            'Isoelectric Point': analysed_seq.isoelectric_point(),
            'Instability Index': analysed_seq.instability_index(),
            'Aromaticity': analysed_seq.aromaticity(),
            'Gravy': analysed_seq.gravy()
        }
    except Exception as e:
        # Fallback for ambiguous sequences or errors
        props = {
            'Molecular Weight': 0,
            'Isoelectric Point': 0,
            'Instability Index': 0,
            'Aromaticity': 0,
            'Gravy': 0,
            'Error': str(e)
        }
        
    return props
