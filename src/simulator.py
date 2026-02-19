from Bio.Seq import Seq
from analysis_advanced import translate_dna
from analysis_protein import calculate_protein_properties

class MutationSimulator:
    def __init__(self, original_dna):
        self.original_dna = str(original_dna).upper()
        self.mutated_dna = str(original_dna).upper()
        
    def apply_mutation(self, position_1based, new_base):
        """
        Apply a point mutation at the given 1-based position.
        """
        idx = position_1based - 1
        if 0 <= idx < len(self.original_dna):
            chars = list(self.original_dna)
            chars[idx] = new_base.upper()
            self.mutated_dna = "".join(chars)
            return True
        return False
        
    def reset(self):
        self.mutated_dna = self.original_dna
        
    def compare_protein_properties(self):
        """
        Compare physicochemical properties of original vs mutated protein.
        """
        prot_orig = translate_dna(self.original_dna)
        prot_mut = translate_dna(self.mutated_dna)
        
        props_orig = calculate_protein_properties(prot_orig)
        props_mut = calculate_protein_properties(prot_mut)
        
        # Calculate diffs
        diffs = {}
        for key in props_orig:
            if isinstance(props_orig[key], (int, float)) and key != 'Error':
                diffs[key] = {
                    'Original': props_orig[key],
                    'Mutated': props_mut[key],
                    'Change': props_mut[key] - props_orig[key]
                }
                
        return {
            'Protein_Original': prot_orig,
            'Protein_Mutated': prot_mut,
            'Properties_Comparison': diffs,
            'Is_Silent': prot_orig == prot_mut,
            'Stop_Codon_Introduced': '*' in prot_mut and '*' not in prot_orig
        }
