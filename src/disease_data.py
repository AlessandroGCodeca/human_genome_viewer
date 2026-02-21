import requests
from typing import List, Dict, Any

def get_disease_variants(gene_id: str) -> List[Dict[str, Any]]:
    """
    Returns a list of disease variants by querying MyVariant.info API (ClinVar data).
    Searches for variants associated with the given Gene Symbol or ID.
    """
    # Clean up the gene_id in case it's an accession
    # MyVariant works best with actual gene symbols (e.g. 'BRCA1', 'INS')
    raw_query = gene_id.split('.')[0] # Remove version if present
    
    # We search the clinvar database for the gene
    url = f"https://myvariant.info/v1/query?q=clinvar.gene.symbol:{raw_query} OR clinvar.gene.id:{raw_query}&fields=clinvar&size=10"
    
    try:
        response = requests.get(url, timeout=5)
        response.raise_for_status()
        data = response.json()
        
        variants = []
        hits = data.get('hits', [])
        
        for hit in hits:
            clinvar = hit.get('clinvar', {})
            if not clinvar:
                continue
                
            # Extract clinical significance
            rcv = clinvar.get('rcv', [])
            sig = "Unknown"
            disease = "Not specified"
            
            if isinstance(rcv, list) and len(rcv) > 0:
                sig = rcv[0].get('clinical_significance', 'Unknown')
                conditions = rcv[0].get('conditions', {})
                if isinstance(conditions, dict):
                    disease = conditions.get('name', 'Not specified')
                elif isinstance(conditions, list) and len(conditions) > 0:
                    disease = conditions[0].get('name', 'Not specified')
            elif isinstance(rcv, dict):
                sig = rcv.get('clinical_significance', 'Unknown')
                conditions = rcv.get('conditions', {})
                if isinstance(conditions, dict):
                    disease = conditions.get('name', 'Not specified')
                    
            # Extract variant representation
            hgvs = clinvar.get('hgvs', {})
            coding = hgvs.get('coding', '')
            genomic = hgvs.get('genomic', '')
            
            desc = coding if coding else genomic
            if not desc:
                desc = f"Variant ID: {clinvar.get('variant_id', 'Unknown')}"
                
            variant_obj = {
                'disease': disease,
                'sig': sig,
                'desc': desc,
                # Real allele freqs are deep in gnomad, providing a placeholder for UI compatibility
                'allele_freq': {'EUR': 'N/A', 'EAS': 'N/A', 'AFR': 'N/A', 'AMR': 'N/A', 'SAS': 'N/A'},
                'freq_desc': 'Frequency data requires gnomAD integration.',
                'pos': clinvar.get('variant_id', 'N/A'),
                'ref': '-',
                'alt': '-'
            }
            variants.append(variant_obj)
            
        return variants
        
    except requests.RequestException as e:
        print(f"Error fetching from MyVariant.info: {e}")
        
    return []
