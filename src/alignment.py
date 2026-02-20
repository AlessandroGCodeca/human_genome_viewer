from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
import matplotlib.pyplot as plt
import io

def format_alignment(alignments):
    """
    Format the best alignment for display.
    """
    if not alignments:
        return "No alignment found."
        
    best_alignment = alignments[0]
    return str(best_alignment)

def perform_pairwise_alignment(seq1, seq2, match_score=1.0, mismatch_score=-1.0, gap_score=-0.5):
    """
    Perform pairwise alignment between two sequences.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = gap_score
    aligner.extend_gap_score = -0.1
    
    alignments = aligner.align(seq1, seq2)
    
    return {
        'score': alignments[0].score if alignments else 0,
        'alignment_str': format_alignment(alignments),
        'count': len(alignments)
    }

def build_phylogenetic_tree(fasta_dict):
    """
    Build a phylogenetic tree from a dictionary of sequences using Pairwise distances.
    
    Args:
        fasta_dict (dict): Dictionary with sequence names as keys and DNA sequences as values.
    
    Returns:
        tuple: (matplotlib Figure, list of distance matrix rows for display)
    """
    names = list(fasta_dict.keys())
    sequences = list(fasta_dict.values())
    n = len(names)
    
    # Initialize lower-triangular distance matrix
    matrix = []
    
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    
    # Calculate pairwise distances
    # Distance = 1.0 - (score / max_score)
    for i in range(n):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                seq1 = sequences[i]
                seq2 = sequences[j]
                
                # Perform alignment
                alignments = aligner.align(seq1, seq2)
                score = alignments[0].score if alignments else 0
                
                # Normalize score to distance (0 to 1)
                # Max possible score is the match_score * length of the shorter sequence
                max_len = max(len(seq1), len(seq2))
                max_score = max_len * aligner.match_score
                distance = 1.0 - (score / max_score) if max_score > 0 else 1.0
                
                # Ensure distance is positive and bounds
                distance = max(0.0, min(1.0, distance))
                row.append(distance)
        matrix.append(row)
    
    # Build distance matrix object
    dm = DistanceMatrix(names, matrix)
    
    # Build Tree using UPGMA
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    
    # Render Tree with Matplotlib
    fig, ax = plt.subplots(figsize=(8, 5))
    Phylo.draw(tree, axes=ax, do_show=False)
    
    # Format matrix for display
    matrix_data = []
    for i, row in enumerate(matrix):
        row_dict = {"Sequence": names[i]}
        for j, dist in enumerate(row):
            row_dict[names[j]] = round(dist, 4)
        matrix_data.append(row_dict)
        
    return fig, matrix_data
