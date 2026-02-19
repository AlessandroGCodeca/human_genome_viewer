import pandas as pd
import os

data_dir = 'data'
files = [f for f in os.listdir(data_dir) if f.endswith('.csv')]

for file in files:
    print(f"--- {file} ---")
    try:
        df = pd.read_csv(os.path.join(data_dir, file))
        print("Columns:", df.columns.tolist())
        print("Shape:", df.shape)
        print("First 2 rows:")
        print(df.head(2))
        
        # Check for sequence-like columns
        for col in df.columns:
            if df[col].dtype == object:
                sample = str(df[col].iloc[0])
                if len(sample) > 20 and all(c in 'ATCGN' for c in sample.upper()):
                    print(f"POTENTIAL SEQUENCE FOUND IN COLUMN: {col}")
    except Exception as e:
        print(f"Error reading {file}: {e}")
    print("\n")
