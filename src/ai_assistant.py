import time

def get_ai_response(gene_id: str, query: str) -> str:
    """
    Mock AI Assistant to explain gene function.
    In a real application, this would call an LLM API.
    """
    query_lower = query.lower()
    
    # Simulate API latency
    time.sleep(1.0)
    
    # Mock responses based on gene IDs
    if 'NM_000014' in gene_id:  # Insulin
        if 'liver' in query_lower:
            return "Insulin signals the liver to take up glucose from the blood and store it as glycogen. It also inhibits the liver from producing new glucose (gluconeogenesis)."
        elif 'function' in query_lower or 'do' in query_lower:
            return "Insulin is a hormone produced by the pancreas. It acts like a key, allowing your body to use sugar (glucose) from carbohydrates in the food that you eat for energy or to store glucose for future use. Insulin helps keeps your blood sugar level from getting too high (hyperglycemia) or too low (hypoglycemia)."
        else:
            return "I am a mock AI. For the Insulin (INS) gene, I can answer questions about its basic function or what it does in the liver."
            
    elif 'NM_007294' in gene_id:  # BRCA1
        if 'cancer' in query_lower:
            return "Mutations in the BRCA1 gene can prevent it from properly repairing damaged DNA, which increases the likelihood that cells will develop additional genetic alterations that can lead to cancer, particularly breast and ovarian cancer."
        elif 'function' in query_lower or 'do' in query_lower:
            return "BRCA1 (Breast Cancer Type 1 susceptibility protein) is a tumor suppressor gene. Its primary role is to repair damaged DNA, playing a crucial role in maintaining the stability of the cell's genetic information."
        else:
            return "I am a mock AI. For the BRCA1 gene, I can answer questions about its role as a tumor suppressor or its link to cancer."
            
    elif 'NM_000518' in gene_id:  # HBB
        if 'sickle' in query_lower or 'disease' in query_lower:
            return "A specific mutation in the HBB gene (Glu6Val) causes red blood cells to form an abnormal sickle or crescent shape. These cells can block blood flow, leading to the symptoms of sickle cell disease."
        elif 'function' in query_lower or 'do' in query_lower:
            return "The HBB gene provides instructions for making a protein called beta-globin. Beta-globin is a subunit of a larger protein called hemoglobin, which is located inside red blood cells. Hemoglobin binds to oxygen in the lungs and carries it to tissues throughout the body."
        else:
            return "I am a mock AI. For the HBB gene, I can answer questions about its normal function in hemoglobin or its relation to sickle cell disease."
            
    else:
        return f"I don't have specific mock information for the gene ID '{gene_id}'. Try asking about Insulin (NM_000014), BRCA1 (NM_007294), or Hemoglobin Beta (NM_000518)."
