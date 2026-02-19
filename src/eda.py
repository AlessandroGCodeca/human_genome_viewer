import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set style
sns.set(style="whitegrid")

DATA_DIR = 'data'
PLOTS_DIR = 'eda_plots'
os.makedirs(PLOTS_DIR, exist_ok=True)

def load_data(filename):
    path = os.path.join(DATA_DIR, filename)
    if os.path.exists(path):
        return pd.read_csv(path)
    return None

def analyze_genomic_summary():
    df = load_data('GRCh38_latest_genomic_summary.csv')
    if df is not None:
        print(f"Genomic Summary shape: {df.shape}")
        
        # Plot sequence length distribution (log scale because huge range)
        plt.figure(figsize=(10, 6))
        sns.histplot(df['Sequence_len'], bins=50, log_scale=True)
        plt.title('Distribution of Genomic Sequence Lengths')
        plt.xlabel('Length (log scale)')
        plt.savefig(os.path.join(PLOTS_DIR, 'genomic_len_dist.png'))
        plt.close()

def analyze_rna_summary():
    df = load_data('GRCh38_latest_rna_summary.csv')
    if df is not None:
        print(f"RNA Summary shape: {df.shape}")
        
        # Plot sequence length distribution
        plt.figure(figsize=(10, 6))
        sns.histplot(df['Sequence_len'], bins=50)
        plt.title('Distribution of RNA Sequence Lengths')
        plt.xlabel('Length')
        plt.savefig(os.path.join(PLOTS_DIR, 'rna_len_dist.png'))
        plt.close()

def analyze_protein_summary():
    df = load_data('GRCh38_latest_protein_symmery.csv') # Note: typo in filename from dataset
    if df is not None:
        print(f"Protein Summary shape: {df.shape}")
        
        # Plot sequence length distribution
        plt.figure(figsize=(10, 6))
        sns.histplot(df['Sequence_len'], bins=50)
        plt.title('Distribution of Protein Sequence Lengths')
        plt.xlabel('Length')
        plt.savefig(os.path.join(PLOTS_DIR, 'protein_len_dist.png'))
        plt.close()

if __name__ == "__main__":
    analyze_genomic_summary()
    analyze_rna_summary()
    analyze_protein_summary()
    print("EDA plots generated in eda_plots/")
