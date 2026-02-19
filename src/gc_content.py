def calculate_gc_content(sequence):
    """
    Calculate the GC content of a DNA sequence.
    
    Args:
        sequence (str or Seq): The DNA sequence.
        
    Returns:
        float: The GC content percentage (0-100).
    """
    if not sequence:
        return 0.0
    
    seq_upper = str(sequence).upper()
    g_count = seq_upper.count('G')
    c_count = seq_upper.count('C')
    total_len = len(seq_upper)
    
    if total_len == 0:
        return 0.0
        
    return ((g_count + c_count) / total_len) * 100.0

if __name__ == "__main__":
    # Test
    test_seq = "ATGCATGC" # 50% GC
    print(f"Sequence: {test_seq}, GC Content: {calculate_gc_content(test_seq)}%")
    
    test_seq2 = "AAAA" # 0% GC
    print(f"Sequence: {test_seq2}, GC Content: {calculate_gc_content(test_seq2)}%")
    
    test_seq3 = "GGCC" # 100% GC
    print(f"Sequence: {test_seq3}, GC Content: {calculate_gc_content(test_seq3)}%")
