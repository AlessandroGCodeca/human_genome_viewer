from src.sequence_fetcher import SequenceFetcher
from src.gc_content import calculate_gc_content
from src.kmer_frequency import calculate_kmer_frequency, get_most_common_kmers

def test_integration():
    print("--- Testing Integration of Sequence Analysis Modules ---")
    
    # 1. Fetch Sequence
    fetcher = SequenceFetcher()
    # Using a short gene for speed, e.g., Human Insulin (INS) gene variant, 
    # but using the one from earlier test (NM_000014.6) is safer as we know it works.
    target_id = "NM_000014.6"
    print(f"Fetching sequence for {target_id}...")
    seq_record = fetcher.fetch_sequence(target_id)
    
    if not seq_record:
        print("FAILED: Could not fetch sequence.")
        return

    sequence = seq_record.seq
    print(f"Successfully fetched sequence of length {len(sequence)}")
    
    # 2. GC Content
    gc = calculate_gc_content(sequence)
    print(f"GC Content: {gc:.2f}%")
    
    # 3. K-mer Frequency
    k = 3
    kmers = calculate_kmer_frequency(sequence, k)
    top_kmers = get_most_common_kmers(kmers, 5)
    print(f"Top 5 {k}-mers:")
    for kmer, count in top_kmers:
        print(f"  {kmer}: {count}")

    print("\n--- Integration Test Complete ---")

if __name__ == "__main__":
    test_integration()
