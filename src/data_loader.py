import pandas as pd
import os

DEFAULT_DATA_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'GRCh38_latest_rna_summary.csv')

def load_gene_metadata(file_path=DEFAULT_DATA_PATH):
    """
    Load the RNA summary CSV file into a pandas DataFrame.
    
    Args:
        file_path (str): Path to the CSV file.
        
    Returns:
        pd.DataFrame: The loaded DataFrame, or None if file not found.
    """
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return None
        
    try:
        # Load only necessary columns to save memory if needed, but current dataset is small enough
        df = pd.read_csv(file_path)
        # Ensure description is string for searching
        df['Description'] = df['Description'].fillna('').astype(str)
        return df
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

def search_genes(df, query, limit=50):
    """
    Search for genes in the DataFrame based on the query.
    
    Args:
        df (pd.DataFrame): The gene metadata DataFrame.
        query (str): The search query string.
        limit (int): Maximum complexity of results to return.
        
    Returns:
        pd.DataFrame: Filtered DataFrame with search results.
    """
    if df is None or df.empty:
        return pd.DataFrame()
    
    if not query:
        return df.head(limit)
    
    # Case-insensitive search in Description or ID
    mask = (
        df['Description'].str.contains(query, case=False, na=False) | 
        df['ID'].str.contains(query, case=False, na=False)
    )
    
    results = df[mask]
    return results.head(limit)
