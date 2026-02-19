from Bio.Seq import Seq
from collections import Counter
import math
import pandas as pd

def translate_dna(sequence):
    """
    Translate a DNA sequence to a Protein sequence.
    
    Args:
        sequence (str or Seq): The DNA sequence.
        
    Returns:
        str: The protein sequence.
    """
    if not sequence:
        return ""
    
    seq_obj = Seq(str(sequence))
    # Using 'to_stop=False' to show full translation including stop codons as '*'
    # If len(seq) is not multiple of 3, Biopython might warn/error, so let's truncate
    trim_len = len(seq_obj) - (len(seq_obj) % 3)
    return str(seq_obj[:trim_len].translate())

def calculate_shannon_entropy(sequence):
    """
    Calculate the Shannon Entropy of a DNA sequence.
    H = -sum(p_i * log2(p_i))
    
    Args:
        sequence (str): The DNA sequence.
        
    Returns:
        float: The entropy value.
    """
    if not sequence:
        return 0.0
        
    seq = str(sequence).upper()
    algo_len = len(seq)
    bases = Counter(seq)
    
    entropy = 0
    for count in bases.values():
        p_i = count / algo_len
        entropy -= p_i * math.log2(p_i)
        
    return entropy

def rolling_entropy(sequence, window_size=100, step=20):
    """
    Calculate entropy over a rolling window.
    
    Args:
        sequence (str): DNA sequence.
        window_size (int): Size of the window.
        step (int): Step size for sliding.
        
    Returns:
        pd.DataFrame: DataFrame with 'Position' and 'Entropy'.
    """
    if not sequence or len(sequence) < window_size:
        return pd.DataFrame()
        
    seq = str(sequence).upper()
    data = []
    
    for i in range(0, len(seq) - window_size + 1, step):
        window_seq = seq[i:i+window_size]
        ent = calculate_shannon_entropy(window_seq)
        data.append({'Position': i, 'Entropy': ent})
        
    return pd.DataFrame(data)

def find_motif(sequence, motif):
    """
    Find all occurrences of a motif in a sequence.
    
    Args:
        sequence (str): DNA sequence.
        motif (str): The pattern to search for (e.g., 'GAATTC').
        
    Returns:
        list: List of starting indices (0-based) where motif is found.
    """
    if not sequence or not motif:
        return []
    
    seq = str(sequence).upper()
    motif = str(motif).upper()
    
    indices = []
    pos = seq.find(motif)
    while pos != -1:
        indices.append(pos)
        pos = seq.find(motif, pos + 1)
        
    return indices
