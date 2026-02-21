import polars as pl
from src.analysis_dna import calculate_gc_skew
from src.analysis_advanced import rolling_entropy

seq = "ATGCGTA"*1000
skew = calculate_gc_skew(seq, window_size=100, step_size=20)
print("Skew passed", len(skew))

entropy = rolling_entropy(seq, window_size=100, step=20)
print("Entropy passed", entropy.shape)
