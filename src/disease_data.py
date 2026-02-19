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
            {'pos': 140, 'ref': 'A', 'alt': 'G', 'disease': 'Diabetes Mellitus', 'sig': 'Pathogenic', 'desc': 'Promoter variant associated with susceptibility.'},
            {'pos': 215, 'ref': 'C', 'alt': 'T', 'disease': 'Hyperproinsulinemia', 'sig': 'Likely Pathogenic', 'desc': 'Missence variant affecting proinsulin processing.'},
            {'pos': 300, 'ref': 'G', 'alt': 'A', 'disease': 'Neonatal Diabetes', 'sig': 'Pathogenic', 'desc': 'Critical region mutation.'}
        ],
        # BRCA1
        'NM_007294': [
            {'pos': 1500, 'ref': 'AG', 'alt': '-', 'disease': 'Breast Cancer Susceptibility', 'sig': 'Pathogenic', 'desc': 'Frameshift deletion.'},
            {'pos': 2300, 'ref': 'T', 'alt': 'G', 'disease': 'Ovarian Cancer Risk', 'sig': 'Pathogenic', 'desc': 'Nonsense mutation.'}
        ],
        # HBB (Sickle Cell)
        'NM_000518': [
            {'pos': 20, 'ref': 'A', 'alt': 'T', 'disease': 'Sickle Cell Anemia', 'sig': 'Pathogenic', 'desc': 'Classic Glu6Val mutation site (approximate loc).'},
            {'pos': 21, 'ref': 'G', 'alt': 'C', 'disease': 'Beta Thalassemia', 'sig': 'Pathogenic', 'desc': 'Splice site variant.'}
        ]
    }
    
    # Simple partial match logic for accession IDs
    for key in db:
        if key in gene_id:
            return db[key]
            
    return []
