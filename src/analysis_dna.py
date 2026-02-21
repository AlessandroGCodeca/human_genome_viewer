import polars as pl
from collections import Counter

def calculate_gc_skew(sequence, window_size=100, step_size=20):
    """
    Calculate GC Skew across a sequence with a sliding window.
    GC Skew = (G - C) / (G + C)
    Optimized with Polars for large scale genomic sequences.
    """
    seq_str = str(sequence).upper()
    
    # Handle short sequences
    if len(seq_str) < window_size:
        g = seq_str.count('G')
        c = seq_str.count('C')
        skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
        return pl.DataFrame({'Position': [0], 'GC_Skew': [skew]}).to_dicts()

    # Create a Polars DataFrame from the sequence
    df = pl.DataFrame({"base": list(seq_str)})
    
    # Add row index for dynamic grouping
    df = df.with_row_index("idx")
    
    # Create indicator columns
    df = df.with_columns(
        (pl.col("base") == "G").cast(pl.Int32).alias("is_g"),
        (pl.col("base") == "C").cast(pl.Int32).alias("is_c")
    )
    
    # Perform rolling sum with step
    rolling_df = df.group_by_dynamic(
        "idx",
        every=f"{step_size}i",
        period=f"{window_size}i",
        include_boundaries=True
    ).agg(
        pl.col("is_g").sum().alias("g_count"),
        pl.col("is_c").sum().alias("c_count")
    )
    
    # Filter out incomplete windows at the end if necessary, though period matches bounds.
    # Calculate Skew
    res_df = rolling_df.with_columns(
        pl.when(pl.col("g_count") + pl.col("c_count") == 0)
        .then(0.0)
        .otherwise((pl.col("g_count") - pl.col("c_count")) / (pl.col("g_count") + pl.col("c_count")))
        .alias("GC_Skew")
    ).select(
        pl.col("idx").alias("Position"),
        pl.col("GC_Skew")
    )
    
    return res_df.to_dicts()

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
