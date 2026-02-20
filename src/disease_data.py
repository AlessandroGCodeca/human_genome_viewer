def get_disease_variants(gene_id):
    """
    Returns a list of mock disease variants for demonstration purposes.
    In a real app, this would query ClinVar or similar API.
    """
    
    # Mock Database
    # Key: Gene Accession ID prefix (or specific ID)
    # Value: List of variants
    
    # NM_000014 is Insulin
    # NM_007294 is BRCA1
    # NM_000518 is HBB (Hemoglobin Beta)
    
    db = {
        # Insulin (INS)
        'NM_000014': [
            {'pos': 140, 'ref': 'A', 'alt': 'G', 'disease': 'Diabetes Mellitus', 'sig': 'Pathogenic', 'desc': 'Promoter variant associated with susceptibility.', 'allele_freq': {'EUR': 0.05, 'EAS': 0.01, 'AFR': 0.02, 'AMR': 0.04, 'SAS': 0.03}, 'freq_desc': 'Common in Europe but rare in Asia.'},
            {'pos': 215, 'ref': 'C', 'alt': 'T', 'disease': 'Hyperproinsulinemia', 'sig': 'Likely Pathogenic', 'desc': 'Missence variant affecting proinsulin processing.', 'allele_freq': {'EUR': 0.001, 'EAS': 0.005, 'AFR': 0.001, 'AMR': 0.002, 'SAS': 0.001}, 'freq_desc': 'Very rare across all populations.'},
            {'pos': 300, 'ref': 'G', 'alt': 'A', 'disease': 'Neonatal Diabetes', 'sig': 'Pathogenic', 'desc': 'Critical region mutation.', 'allele_freq': {'EUR': 0.0001, 'EAS': 0.0001, 'AFR': 0.0001, 'AMR': 0.0001, 'SAS': 0.0001}, 'freq_desc': 'Extremely rare de novo mutation.'}
        ],
        # BRCA1
        'NM_007294': [
            {'pos': 1500, 'ref': 'AG', 'alt': '-', 'disease': 'Breast Cancer Susceptibility', 'sig': 'Pathogenic', 'desc': 'Frameshift deletion.', 'allele_freq': {'EUR': 0.01, 'EAS': 0.001, 'AFR': 0.005, 'AMR': 0.008, 'SAS': 0.002}, 'freq_desc': 'Higher prevalence in individuals of European descent.'},
            {'pos': 2300, 'ref': 'T', 'alt': 'G', 'disease': 'Ovarian Cancer Risk', 'sig': 'Pathogenic', 'desc': 'Nonsense mutation.', 'allele_freq': {'EUR': 0.005, 'EAS': 0.002, 'AFR': 0.001, 'AMR': 0.003, 'SAS': 0.001}, 'freq_desc': 'Rare across all populations.'}
        ],
        # HBB (Sickle Cell)
        'NM_000518': [
            {'pos': 20, 'ref': 'A', 'alt': 'T', 'disease': 'Sickle Cell Anemia', 'sig': 'Pathogenic', 'desc': 'Classic Glu6Val mutation site (approximate loc).', 'allele_freq': {'EUR': 0.001, 'EAS': 0.001, 'AFR': 0.15, 'AMR': 0.05, 'SAS': 0.02}, 'freq_desc': 'Highly prevalent in African populations due to malaria protection.'},
            {'pos': 21, 'ref': 'G', 'alt': 'C', 'disease': 'Beta Thalassemia', 'sig': 'Pathogenic', 'desc': 'Splice site variant.', 'allele_freq': {'EUR': 0.02, 'EAS': 0.01, 'AFR': 0.03, 'AMR': 0.01, 'SAS': 0.08}, 'freq_desc': 'Common in Mediterranean and South Asian populations.'}
        ]
    }
    
    # Simple partial match logic for accession IDs
    for key in db:
        if key in gene_id:
            return db[key]
            
    return []
