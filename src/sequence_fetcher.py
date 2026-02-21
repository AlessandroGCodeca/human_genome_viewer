import os
from Bio import Entrez, SeqIO
from urllib.error import HTTPError
import time
from typing import List, Dict, Optional, Tuple, Any
from tenacity import retry, wait_exponential, stop_after_attempt, retry_if_exception_type
from concurrent.futures import ThreadPoolExecutor, as_completed

class SequenceFetcher:
    def __init__(self, email: str = "tool_user@example.com", retmax: int = 1):
        """
        Initialize the SequenceFetcher.
        
        Args:
            email (str): Email address is required by NCBI Entrez.
            retmax (int): Maximum number of records to retrieve per request.
        """
        Entrez.email = email
        self.retmax = retmax

    @retry(wait=wait_exponential(multiplier=1, min=2, max=10), stop=stop_after_attempt(3))
    def _fetch_sequence_with_retry(self, accession_id: str, db: str) -> Any:
        handle = Entrez.efetch(db=db, id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record

    def fetch_sequence(self, accession_id: str, db: str = "nucleotide") -> Optional[Any]:
        """
        Fetch a sequence from NCBI with exponential backoff.

        Args:
            accession_id (str): The accession ID (e.g., 'NM_000014.6').
            db (str): The database to search in (default: 'nucleotide').

        Returns:
            SeqRecord: A Biopython SeqRecord object, or None if failed.
        """
        try:
            print(f"Fetching {accession_id} from {db}...")
            return self._fetch_sequence_with_retry(accession_id, db)
        except Exception as e:
            print(f"Error fetching {accession_id}: {e}")
            return None
        finally:
            time.sleep(0.34) # Be nice to NCBI servers

    @retry(wait=wait_exponential(multiplier=1, min=2, max=10), stop=stop_after_attempt(3))
    def _esearch_with_retry(self, db: str, term: str, retmax: int) -> dict:
        handle = Entrez.esearch(db=db, term=term, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        return record

    @retry(wait=wait_exponential(multiplier=1, min=2, max=10), stop=stop_after_attempt(3))
    def _esummary_with_retry(self, db: str, id_str: str) -> List[dict]:
        summary_handle = Entrez.esummary(db=db, id=id_str)
        summaries = Entrez.read(summary_handle)
        summary_handle.close()
        return summaries

    def search_gene_by_name(self, gene_name: str, organism: str = "Homo sapiens", limit: int = 10) -> List[Dict[str, str]]:
        """
        Search for a gene by name in the nucleotide database using NCBI E-utilities.
        
        Args:
            gene_name (str): The name of the gene (e.g., 'TP53').
            organism (str): The organism to filter by.
            limit (int): Max results to return.
            
        Returns:
            list: A list of dictionaries containing 'ID' and 'Description'.
        """
        try:
            term = f"({gene_name}[Gene]) AND {organism}[Organism] AND mRNA[Filter] AND refseq[Filter]"
            record = self._esearch_with_retry(db="nucleotide", term=term, retmax=limit)
            time.sleep(0.34)
            
            id_list = record.get("IdList", [])
            if not id_list:
                return []
                
            summaries = self._esummary_with_retry(db="nucleotide", id_str=",".join(id_list))
            time.sleep(0.34)
            
            results = []
            for summary in summaries:
                acc_version = summary.get('AccessionVersion', '')
                title = summary.get('Title', '')
                if acc_version:
                    results.append({
                        'ID': acc_version,
                        'Description': title
                    })
            return results
        except Exception as e:
            print(f"Error searching gene {gene_name}: {e}")
            return []

    def fetch_gene_location(self, accession_id: str) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
        """
        Fetches the chromosome and cytoband location for a given accession ID.
        This requires a two-step process: find the Gene ID in the nucleotide DB, then
        query the gene DB for the exact MapLocation.
        """
        try:
            # 1. First search the gene database using the accession ID
            search_results = self._esearch_with_retry(db="gene", term=accession_id, retmax=1)
            time.sleep(0.34)
            
            if not search_results.get("IdList"):
                return None, None, None, None
                
            gene_id = search_results["IdList"][0]
            
            # 2. Fetch the summary for that specific Gene ID
            summary_results = self._esummary_with_retry(db="gene", id_str=gene_id)
            time.sleep(0.34)
            
            if summary_results and len(summary_results["DocumentSummarySet"]["DocumentSummary"]) > 0:
                doc = summary_results["DocumentSummarySet"]["DocumentSummary"][0]
                chrom = doc.get("Chromosome", "")
                map_loc = doc.get("MapLocation", "")
                
                start = None
                stop = None
                genomic_info = doc.get("GenomicInfo", [])
                if genomic_info and len(genomic_info) > 0:
                    info = genomic_info[0]
                    start = info.get("ChrStart")
                    stop = info.get("ChrStop")
                    
                if chrom and map_loc:
                    return str(chrom), str(map_loc), str(start) if start is not None else None, str(stop) if stop is not None else None
            
            return None, None, None, None
            
        except Exception as e:
            print(f"Error fetching location for {accession_id}: {e}")
            return None, None, None, None

    def fetch_multiple_sequences(self, accession_ids: List[str], db: str = "nucleotide") -> Dict[str, Any]:
        """
        Fetch multiple sequences concurrently.
        """
        results = {}
        with ThreadPoolExecutor(max_workers=3) as executor:
            future_to_id = {executor.submit(self.fetch_sequence, acc_id, db): acc_id for acc_id in accession_ids}
            for future in as_completed(future_to_id):
                acc_id = future_to_id[future]
                try:
                    record = future.result()
                    results[acc_id] = record
                except Exception as exc:
                    print(f'{acc_id} generated an exception: {exc}')
                    results[acc_id] = None
        return results

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
