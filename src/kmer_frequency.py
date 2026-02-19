from collections import Counter

def calculate_kmer_frequency(sequence, k=3):
    """
    Calculate the frequency of K-mers in a DNA sequence.
    
    Args:
        sequence (str or Seq): The DNA sequence.
        k (int): Length of the K-mer.
        
    Returns:
        dict: Dictionary where keys are K-mers and values are their counts.
    """
    if not sequence or k <= 0:
        return {}
    
    seq_upper = str(sequence).upper()
    kmers = [seq_upper[i:i+k] for i in range(len(seq_upper) - k + 1)]
    
    return Counter(kmers)

def get_most_common_kmers(frequency_dict, n=10):
    """
    Get the n most common K-mers from a frequency dictionary.
    
    Args:
        frequency_dict (dict): Output from calculate_kmer_frequency.
        n (int): Number of top K-mers to return.
        
    Returns:
        list: List of (kmer, count) tuples.
    """
    if not frequency_dict:
        return []
    
    # Counter object handles most_common directly if passed a Counter, 
    # but here we accept a generic dict
    if isinstance(frequency_dict, Counter):
        return frequency_dict.most_common(n)
    else:
        sorted_kmers = sorted(frequency_dict.items(), key=lambda item: item[1], reverse=True)
        return sorted_kmers[:n]

if __name__ == "__main__":
    # Test
    test_seq = "ATGCATGCATGC" # "ATG", "TGC", "GCA", "CAT"
    k = 3
    freq = calculate_kmer_frequency(test_seq, k)
    print(f"Sequence: {test_seq}, K={k}")
    print("Frequencies:", freq)
    
    print("Most common:", get_most_common_kmers(freq, 2))
