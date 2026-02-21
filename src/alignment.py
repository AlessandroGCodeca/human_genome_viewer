from Bio.Align import PairwiseAligner, MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
import plotly.graph_objects as go
import io
import networkx as nx
import subprocess
import tempfile
import os
from typing import Dict, List, Any, Tuple, Optional
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

def format_alignment(alignments: Any) -> str:
    """
    Format the best alignment for display.
    """
    if not alignments:
        return "No alignment found."
        
    best_alignment = alignments[0]
    return str(best_alignment)

def perform_pairwise_alignment(seq1: str, seq2: str, match_score: float = 1.0, mismatch_score: float = -1.0, gap_score: float = -0.5, mode: str = 'global') -> Dict[str, Any]:
    """
    Perform pairwise alignment between two sequences (Global or Local).
    """
    aligner = PairwiseAligner()
    aligner.mode = mode if mode in ['global', 'local'] else 'global'
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = gap_score
    aligner.extend_gap_score = gap_score / 2.0  # Common heuristic
    
    alignments = aligner.align(seq1, seq2)
    
    return {
        'score': alignments[0].score if alignments else 0,
        'alignment_str': format_alignment(alignments),
        'count': len(alignments),
        'mode': aligner.mode
    }

def build_phylogenetic_tree(fasta_dict: Dict[str, str]) -> Tuple[Optional[go.Figure], List[Dict[str, Any]]]:
    """
    Build an interactive phylogenetic tree from a dictionary of sequences using Pairwise distances.
    
    Args:
        fasta_dict (dict): Dictionary with sequence names as keys and DNA sequences as values.
    
    Returns:
        tuple: (Plotly Figure, list of distance matrix rows for display)
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
                
                # Perform alignment and get only the score (avoids path generation OverflowError)
                score = aligner.score(seq1, seq2)
                
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
    
    # Render Tree with Plotly (Convert Biopython tree to NetworkX for coordinates)
    net = Phylo.to_networkx(tree)
    
    # We use roughly a spectral layout or kamada_kawai for unrooted viewing
    pos = nx.spring_layout(net, seed=42) 
    
    edge_x = []
    edge_y = []
    for edge in net.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    node_text = []
    
    for node in net.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        # Display name if it's a leaf node, otherwise a generic junction point
        name = getattr(node, 'name', '')
        node_text.append(f"Clade: {name}" if name else "Junction Node")

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        text=[name if getattr(node, 'name', None) and "Inner" not in name else "" for node in net.nodes()],
        textposition="top center",
        marker=dict(
            showscale=False,
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            line_width=2))

    fig = go.Figure(data=[edge_trace, node_trace],
             layout=go.Layout(
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                annotations=[ dict(
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.005, y=-0.002 ) ],
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
             )
    
    # Format matrix for display
    matrix_data = []
    for i, row in enumerate(matrix):
        row_dict = {"Sequence": names[i]}
        for j, dist in enumerate(row):
            row_dict[names[j]] = round(dist, 4)
        matrix_data.append(row_dict)
        
    return fig, matrix_data

def perform_multiple_sequence_alignment(fasta_dict: Dict[str, str]) -> str:
    """
    Perform True MSA using 'mafft' via subprocess. 
    If mafft is executing, reads its stdout.
    If mafft is missing (e.g., local machine missing apt-get), falls back to the progressive heuristic.
    """
    if len(fasta_dict) < 2:
        return "Need at least 2 sequences for alignment."

    # Try True MSA
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_in:
            for name, seq in fasta_dict.items():
                tmp_in.write(f">{name}\n{seq}\n")
            tmp_in_path = tmp_in.name

        # Call MAFFT (assuming it's in PATH via packages.txt / apt-get)
        result = subprocess.run(
            ['mafft', '--quiet', '--auto', tmp_in_path],
            capture_output=True,
            text=True,
            check=True
        )
        os.remove(tmp_in_path)
        
        # Format the fasta stdout output slightly for UI readability
        lines = result.stdout.split('\n')
        output = f"--- True Multiple Sequence Alignment (MAFFT) ---\n"
        # We don't want to spit out 10,000 characters of raw multiline FASTA in a text block, 
        # so let's summarize or just return the raw alignments block.
        output += "\n".join(lines[:100]) # Cap at 100 lines for visual safety
        if len(lines) > 100:
            output += "\n... [Truncated for display] ..."
        return output

    except (subprocess.CalledProcessError, FileNotFoundError):
        # Fallback to Progressive Heuristic if mafft is not installed
        res = "⚠️ 'mafft' aligner not found in environment. Falling back to Basic Progressive Heuristic...\n\n"
        res += "--- Basic Progressive Alignment Synopsis ---\n"
        
        sequences = list(fasta_dict.values())
        names = list(fasta_dict.keys())
        res += f"Aligned {len(names)} sequences based on pairwise progressive scoring.\n"
        
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        
        # Demo heuristic: Align all against the first sequence
        base_seq = sequences[0]
        res += f"Reference: {names[0]} (Length {len(base_seq)})\n\n"
        
        for i in range(1, len(sequences)):
            alignments = aligner.align(base_seq, sequences[i])
            best = alignments[0]
            res += f"> {names[i]} (Score: {best.score})\n"
            # Just show the first 100 char summarized snippet for readability
            snippet = str(best).split("\n")
            res += f"{snippet[0][:100]}...\n{snippet[1][:100]}...\n{snippet[2][:100]}...\n\n"
            
        return res
