import unittest
import sys
import os

# Add the parent directory to the sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '.')))

from BioML.main import Sequence, Invalid


class TestSequence(unittest.TestCase):

    def test_valid_dna_sequence(self):
        seq = Sequence("ATGC")
        self.assertEqual(seq.sequence, "ATGC")

    def test_invalid_dna_sequence(self):
        with self.assertRaises(Invalid):
            Sequence("ATGX")

    def test_valid_rna_sequence(self):
        seq = Sequence("AUGC", is_rna=True)
        self.assertEqual(seq.sequence, "AUGC")

    def test_invalid_rna_sequence(self):
        with self.assertRaises(Invalid):
            Sequence("AUGX", is_rna=True)

    def test_reverse_dna_sequence(self):
        seq = Sequence("ATGC")
        self.assertEqual(seq.revseq(), "CGTA")

    def test_reverse_rna_sequence(self):
        seq = Sequence("AUGC", is_rna=True)
        self.assertEqual(seq.revseq(), "CGUA")

    def test_complement_dna_sequence(self):
        seq = Sequence("ATGC")
        self.assertEqual(seq.complseq(), "TACG")

    def test_complement_rna_sequence(self):
        seq = Sequence("AUGC", is_rna=True)
        self.assertEqual(seq.complseq(), "UACG")

    def test_reverse_complement_dna_sequence(self):
        seq = Sequence("ATGC")
        self.assertEqual(seq.complrevseq(), "GCAT")

    def test_reverse_complement_rna_sequence(self):
        seq = Sequence("AUGC", is_rna=True)
        self.assertEqual(seq.complrevseq(), "GCAU")

    def test_str_representation(self):
        seq = Sequence("ATGC")
        self.assertEqual(str(seq), "ATGC")
        
    def test_countnucs(self):
        seq = Sequence("AAGGTTCC")
        self.assertEqual(seq.countnucs(), {'A': 2, 'C': 2, 'G': 2, 'T': 2, 'U':0})
        rna_seq = Sequence("AAGGUUCC", is_rna=True)
        self.assertEqual(rna_seq.countnucs(), {'A': 2, 'C': 2, 'G': 2, 'T':0, 'U': 2})
        with self.assertRaises(Invalid):
            invalid_seq = Sequence("AAGGXTCC")
            invalid_seq.countnucs()
        
    def test_gc_content(self):
        seq = Sequence("AAGGTTCC")
        self.assertAlmostEqual(seq.gc_content(), 50.0)
        seq_empty = Sequence("")
        self.assertAlmostEqual(seq_empty.gc_content(), 0.0)
        
    def test_transDNAtoRNA(self):
        seq = Sequence("ATGC")
        rna_seq = seq.transDNAtoRNA()
        self.assertTrue(rna_seq.is_rna)
        self.assertEqual(str(rna_seq), "AUGC")
        with self.assertRaises(ValueError):
            rna_seq.transDNAtoRNA()

    def test_translate_dna_to_protein(self):
        seq = Sequence("ATGCCCTAG")
        self.assertEqual(seq.transDNAtoP(), "MP")

        seq_with_stop = Sequence("ATGTAA")
        self.assertEqual(seq_with_stop.transDNAtoP(), "M")
        
        with self.assertRaises(Invalid):
            rna_seq = Sequence("AUGC", is_rna=True)
            rna_seq.transDNAtoP()

    def test_translate_dna_all_frames(self):
        seq = Sequence("ATGCCCTAG")
        self.assertEqual(seq.transDNAtoP_all(), ["MP", "CP", "AL"])

        with self.assertRaises(Invalid):
            rna_seq = Sequence("AUGC", is_rna=True)
            rna_seq.transDNAtoP_all()
            
    def test_transRNAtoP(self):
        rna_seq = Sequence("AUGCCCUAG", is_rna=True)
        self.assertEqual(rna_seq.transRNAtoP(), "MP")

        rna_seq_with_stop = Sequence("AUGUAA", is_rna=True)
        self.assertEqual(rna_seq_with_stop.transRNAtoP(), "M")

        empty_rna_seq = Sequence("", is_rna=True)
        self.assertEqual(empty_rna_seq.transRNAtoP(), "")   
if __name__ == '__main__':
    unittest.main()