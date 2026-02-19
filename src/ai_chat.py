import random

class AIGeneAssistant:
    def __init__(self):
        self.gene_knowledge = {
            'NM_000014': { # Insulin
                'name': 'Insulin (INS)',
                'function': 'Regulates blood glucose levels by facilitating cellular glucose uptake.',
                'disease': 'Defects are associated with Diabetes Mellitus.',
                'mechanism': 'It is a peptide hormone produced by beta cells of the pancreatic islets.',
            },
            'NM_007294': { # BRCA1
                'name': 'BRCA1 DNA Repair Associated',
                'function': 'A tumor suppressor gene that helps repair damaged DNA and maintain genomic stability.',
                'disease': 'Mutations increase the risk of breast and ovarian cancers.',
                'mechanism': 'It acts as a caretaker of the genome by repairing double-strand breaks.',
            },
            'NM_000518': { # HBB
                'name': 'Hemoglobin Subunit Beta',
                'function': 'Provides instructions for making beta-globin, part of hemoglobin.',
                'disease': 'Mutations cause Sickle Cell Anemia and Beta-Thalassemia.',
                'mechanism': 'Hemoglobin transports oxygen from the lungs to tissues.',
            },
            'NM_000546': { # TP53
                'name': 'Tumor Protein p53',
                'function': 'The "Guardian of the Genome". It regulates cell division and prevents mutations.',
                'disease': 'Associated with Li-Fraumeni Syndrome and many sporadic cancers.',
                'mechanism': 'It activates DNA repair proteins or initiates apoptosis if DNA damage is irreparable.',
            }
        }
        
    def get_response(self, query, accession_id):
        """
        Generate a simulated AI response based on the gene context and user query.
        """
        context = self.gene_knowledge.get(accession_id.split('.')[0], None)
        query_lower = query.lower()
        
        if not context:
            return "I'm still learning about this specific gene. My knowledge base currently covers Insulin, BRCA1, HBB, and TP53."
            
        # Simulated "Intelligence" via keyword matching
        if 'function' in query_lower or 'do' in query_lower or 'what is' in query_lower:
            return f"**Function of {context['name']}**: {context['function']}"
            
        if 'disease' in query_lower or 'cancer' in query_lower or 'risk' in query_lower or 'sickle' in query_lower:
            return f"**Disease Association**: {context['disease']}"
            
        if 'mechanism' in query_lower or 'how' in query_lower or 'work' in query_lower:
            return f"**Mechanism**: {context['mechanism']}"
            
        if 'summary' in query_lower or 'tell me about' in query_lower:
            return (f"**{context['name']}**\n\n"
                    f"• {context['function']}\n"
                    f"• {context['mechanism']}\n"
                    f"• {context['disease']}")
            
        # Fallback response
        return (f"I can tell you about the **function**, **mechanism**, or **disease risks** of {context['name']}. "
                "What would you like to know?")
