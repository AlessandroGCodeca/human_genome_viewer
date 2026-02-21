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

    def search_gene_by_name(self, gene_name, organism="Homo sapiens", limit=10):
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
            handle = Entrez.esearch(db="nucleotide", term=term, retmax=limit)
            record = Entrez.read(handle)
            handle.close()
            time.sleep(0.34)
            
            id_list = record.get("IdList", [])
            if not id_list:
                return []
                
            summary_handle = Entrez.esummary(db="nucleotide", id=",".join(id_list))
            summaries = Entrez.read(summary_handle)
            summary_handle.close()
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

    def fetch_gene_location(self, accession_id):
        """
        Fetches the chromosome and cytoband location for a given accession ID.
        This requires a two-step process: find the Gene ID in the nucleotide DB, then
        query the gene DB for the exact MapLocation.
        """
        try:
            # 1. First search the gene database using the accession ID
            search_handle = Entrez.esearch(db="gene", term=accession_id)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            time.sleep(0.34)
            
            if not search_results["IdList"]:
                return None, None
                
            gene_id = search_results["IdList"][0]
            
            # 2. Fetch the summary for that specific Gene ID
            summary_handle = Entrez.esummary(db="gene", id=gene_id)
            summary_results = Entrez.read(summary_handle)
            summary_handle.close()
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
                    return str(chrom), str(map_loc), start, stop
            
            return None, None, None, None
            
        except HTTPError as e:
            print(f"HTTP Error fetching location for {accession_id}: {e}")
            return None, None
        except Exception as e:
            print(f"Error fetching location for {accession_id}: {e}")
            return None, None

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
