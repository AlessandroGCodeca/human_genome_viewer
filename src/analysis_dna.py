from Bio.SeqUtils import GC
from collections import Counter

def calculate_gc_skew(sequence, window_size=100, step_size=20):
    """
    Calculate GC Skew across a sequence with a sliding window.
    GC Skew = (G - C) / (G + C)
    """
    skews = []
    positions = []
    
    seq_str = str(sequence).upper()
    
    # Handle short sequences
    if len(seq_str) < window_size:
        g = seq_str.count('G')
        c = seq_str.count('C')
        if g + c == 0:
            skew = 0.0
        else:
            skew = (g - c) / (g + c)
        return [{'Position': 0, 'GC_Skew': skew}]

    for i in range(0, len(seq_str) - window_size + 1, step_size):
        window = seq_str[i:i+window_size]
        g = window.count('G')
        c = window.count('C')
        
        if g + c == 0:
            skew = 0.0
        else:
            skew = (g - c) / (g + c)
            
        skews.append(skew)
        positions.append(i)
        
    return [{'Position': pos, 'GC_Skew': skew} for pos, skew in zip(positions, skews)]

def calculate_codon_usage(sequence):
    """
    Calculate the frequency of each codon in the sequence.
    """
    seq_str = str(sequence).upper()
    # Ensure sequence length is multiple of 3 for simple codon counting
    # In real app, might want to find ORF, but distinct codons are fine for stats
    
    codons = [seq_str[i:i+3] for i in range(0, len(seq_str), 3) if len(seq_str[i:i+3]) == 3]
    codon_counts = Counter(codons)
    total_codons = len(codons)
    
    usage = []
    for codon, count in codon_counts.items():
        usage.append({
            'Codon': codon,
            'Count': count,
            'Frequency': count / total_codons if total_codons > 0 else 0
        })
        
    # Sort by Amino Acid if we had a table, but for now sort by Codon
    usage.sort(key=lambda x: x['Codon'])
    
    return usage
