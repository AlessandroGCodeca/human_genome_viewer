from src.analysis_advanced import translate_dna, calculate_shannon_entropy, find_motif

def test_features():
    print("Testing Advanced Analysis Features...")
    
    # Translation
    dna = "ATGGCC" # Met-Ala
    prot = translate_dna(dna)
    print(f"DNA: {dna} -> Protein: {prot}")
    assert prot == "MA" or prot == "M A", "Translation failed"
    
    # Entropy
    # High complexity
    dna_high = "ACGTACGT" # Uniform distribution
    ent_high = calculate_shannon_entropy(dna_high)
    print(f"Entropy of {dna_high}: {ent_high:.4f}")
    
    # Low complexity
    dna_low = "AAAA" 
    ent_low = calculate_shannon_entropy(dna_low)
    print(f"Entropy of {dna_low}: {ent_low:.4f}")
    assert ent_low == 0.0, "Entropy of AAAA should be 0"
    
    # Motif Search
    dna_motif = "GAATTCGAATTC"
    indices = find_motif(dna_motif, "GAATTC")
    print(f"Motif 'GAATTC' found at: {indices}")
    assert indices == [0, 6], "Motif search failed"
    
    print("All tests passed!")

if __name__ == "__main__":
    test_features()
