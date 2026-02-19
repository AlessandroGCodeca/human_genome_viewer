import os
from Bio import Entrez, SeqIO
from urllib.error import HTTPError
import time

class SequenceFetcher:
    def __init__(self, email="tool_user@example.com", retmax=1):
        """
        Initialize the SequenceFetcher.
        
        Args:
            email (str): Email address is required by NCBI Entrez.
            retmax (int): Maximum number of records to retrieve per request.
        """
        Entrez.email = email
        self.retmax = retmax

    def fetch_sequence(self, accession_id, db="nucleotide"):
        """
        Fetch a sequence from NCBI.

        Args:
            accession_id (str): The accession ID (e.g., 'NM_000014.6').
            db (str): The database to search in (default: 'nucleotide').

        Returns:
            SeqRecord: A Biopython SeqRecord object, or None if failed.
        """
        try:
            print(f"Fetching {accession_id} from {db}...")
            # efetch retrieves the full record
            handle = Entrez.efetch(db=db, id=accession_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            return record
        except HTTPError as e:
            print(f"HTTP Error fetching {accession_id}: {e}")
            return None
        except Exception as e:
            print(f"Error fetching {accession_id}: {e}")
            return None
        finally:
            # Be nice to NCBI servers
            time.sleep(0.34) # limit to 3 requests per second

if __name__ == "__main__":
    # Test
    fetcher = SequenceFetcher()
    # Test with a known ID from the RNA summary (NM_000014.6 was in the inspection output)
    seq = fetcher.fetch_sequence("NM_000014.6")
    if seq:
        print(f"Success! Fetched {seq.id}")
        print(f"Length: {len(seq)}")
        print(f"First 100 bases: {seq.seq[:100]}")
    else:
        print("Failed to fetch sequence.")
