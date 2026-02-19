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
            {'pos': 140, 'ref': 'A', 'alt': 'G', 'disease': 'Diabetes Mellitus', 'sig': 'Pathogenic', 'desc': 'Missence mutation affecting insulin stability.', 'frequencies': {'EUR': 0.05, 'AFR': 0.01, 'EAS': 0.00, 'SAS': 0.02, 'AMR': 0.04}},
            {
            'pos': 210, 
            'ref': 'C', 
            'alt': 'T', 
            'disease': 'Hyperproinsulinemia', 
            'sig': 'Likely Pathogenic',
            'desc': 'Disrupts cleavage site.',
            'frequencies': {'EUR': 0.01, 'AFR': 0.08, 'EAS': 0.02, 'SAS': 0.01, 'AMR': 0.02}
        }
    ],
    'NM_007294': [ # BRCA1
        {
            'pos': 1500, 
            'ref': 'AG', 
            'alt': '-', 
            'disease': 'Breast Cancer Susceptibility', 
            'sig': 'Pathogenic',
            'desc': 'Frameshift deletion leading to premature stop.',
            'frequencies': {'EUR': 0.02, 'AFR': 0.00, 'EAS': 0.00, 'SAS': 0.01, 'AMR': 0.01}
        },
        {
            'pos': 2300, 
            'ref': 'T', 
            'alt': 'G', 
            'disease': 'Ovarian Cancer Risk', 
            'sig': 'Pathogenic',
            'desc': 'Nonsense mutation.',
            'frequencies': {'EUR': 0.03, 'AFR': 0.01, 'EAS': 0.01, 'SAS': 0.04, 'AMR': 0.02}
        }],
        # HBB (Sickle Cell)
        'NM_000518': [ # HBB
        {
            'pos': 60, 
            'ref': 'A', 
            'alt': 'T', 
            'disease': 'Sickle Cell Anemia', 
            'sig': 'Pathogenic',
            'desc': 'Glu6Val substitution (HbS).',
            'frequencies': {'EUR': 0.01, 'AFR': 0.15, 'EAS': 0.00, 'SAS': 0.03, 'AMR': 0.05}
        }
    ],
    'NM_000546': [ # TP53
        {
            'pos': 700, 
            'ref': 'G', 
            'alt': 'A', 
            'disease': 'Li-Fraumeni Syndrome', 
            'sig': 'Pathogenic',
            'desc': 'Loss of function in tumor suppressor.',
            'frequencies': {'EUR': 0.00, 'AFR': 0.00, 'EAS': 0.01, 'SAS': 0.00, 'AMR': 0.00}
        }
    ]}
    
    # Simple partial match logic for accession IDs
    for key in db:
        if key in gene_id:
            return db[key]
            
    return []
