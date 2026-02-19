from Bio.Align import PairwiseAligner

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
