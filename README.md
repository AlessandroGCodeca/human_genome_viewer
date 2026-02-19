# Human Genome Analysis Dashboard

An interactive Streamlit application for analyzing sequence data from the Human Genome Dataset. This tool fetches DNA sequences from NCBI, performs bioinformatics analysis (GC content, K-mer frequency, translation, entropy), and visualizes the results.

## Features

-   **Gene Search:** Search for genes (e.g., "Insulin", "TP53") using local metadata.
-   **Sequence Fetching:** Retrieve sequences directly from NCBI using Accession IDs.
-   **Analysis Tools:**
    -   **GC Content:** Calculate base composition.
    -   **K-mer Frequency:** Identify recurring sequence motifs.
    -   **Sequence Translation:** Convert DNA to Protein sequences.
    -   **Advanced Analysis:** Shannon Entropy (complexity) and Motif Search.
-   **Caching:** Optimized for performance with Streamlit caching.

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/<your-username>/human_genome_viewer.git
    cd human_genome_viewer
    ```

2.  **Set up Virtual Environment:**
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```

3.  **Install Dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

4.  **Download Data:**
    -   Obtain the `aliabedimadiseh/human-genome-dataset` from Kaggle.
    -   Place the summary CSV files (`GRCh38_latest_rna_summary.csv` etc.) into the `data/` directory.

## Usage

Run the Streamlit app:
```bash
streamlit run src/app.py
```

## Exploratory Data Analysis (EDA)

See `eda_report.md` (generated locally) for initial insights into the dataset metadata.

## License

MIT
