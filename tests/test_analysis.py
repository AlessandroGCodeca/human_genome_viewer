import unittest
from src.analysis_dna import calculate_gc_skew, calculate_codon_usage
from src.analysis_protein import calculate_protein_properties

class TestGenomicAnalysis(unittest.TestCase):

    def test_gc_skew(self):
        # Test 1: Simple sequence
        seq = "GGCC" * 10 
        # G=20, C=20 -> Skew = 0
        skew = calculate_gc_skew(seq, window_size=40, step_size=40)
        self.assertEqual(len(skew), 1)
        self.assertAlmostEqual(skew[0]['GC_Skew'], 0.0)

        # Test 2: G rich
        seq_g = "GGGG" * 10
        skew_g = calculate_gc_skew(seq_g, window_size=40)
        self.assertAlmostEqual(skew_g[0]['GC_Skew'], 1.0) # (40-0)/40 = 1

        # Test 3: Short sequence
        res = calculate_gc_skew("A", window_size=10)
        self.assertEqual(res[0]['GC_Skew'], 0.0)

    def test_codon_usage(self):
        # ATG GCC ATG
        seq = "ATGGCCATG" 
        usage = calculate_codon_usage(seq)
        
        # ATG: 2, GCC: 1
        atg = next(item for item in usage if item['Codon'] == 'ATG')
        self.assertEqual(atg['Count'], 2)
        self.assertAlmostEqual(atg['Frequency'], 2/3)
        
        gcc = next(item for item in usage if item['Codon'] == 'GCC')
        self.assertEqual(gcc['Count'], 1)

    def test_protein_properties(self):
        # Example protein sequence
        # M: Met, A: Ala
        seq = "MA" 
        props = calculate_protein_properties(seq)
        
        # Molecular weights: M ~149.2, A ~89.1 (minus water mass in bond)
        # BioPython calc is precise, just check it returns numbers
        self.assertGreater(props['Molecular Weight'], 0)
        self.assertTrue('Isoelectric Point' in props)
        self.assertTrue('Instability Index' in props)
        
        # Empty/Invalid
        props_empty = calculate_protein_properties("")
        self.assertEqual(props_empty['Molecular Weight'], 0)

if __name__ == '__main__':
    unittest.main()
