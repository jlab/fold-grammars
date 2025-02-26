import sys
from os.path import dirname, join
sys.path.append(dirname(__file__) + '/../')

from unittest import TestCase, main

from input import read_fasta, read_CT_file, disentangle_knots

def getFP(filepath):
    return join(dirname(__file__), filepath)

class TestExecute(TestCase):
    def setUp(self):
        self.fp_fasta_mirna = getFP('data/mirnas.fasta')
        self.fp_multiCT = getFP('data/multi.ct')
        self.fp_sars = getFP('data/SARS-CoV-2_Full_Length_Secondary_Structure_Map.ct')

    def tearDown(self):
        pass

    def test_read_fasta(self):
        exp = [
            ('hsa-miR-21-5p::NC_000017.11:59841273-59841295(+)', 'AGCTTATCAGACTGATGTTGAC'),
            ('hsa-let-7i-5p::NC_000012.12:62603691-62603713(+)', 'GAGGTAGTAGTTTGTGCTGTTG'),
            ('hsa-miR-192-5p::NC_000011.10:64891203-64891224(-)', 'TCTGACCTATGAATTGACAGC'),
            ('hsa-miR-151a-3p::NC_000008.11:140732587-140732608(-)', 'ACTAGACtgaagctccttgag'),
            ('hsa-miR-27b-3p::NC_000009.12:95085505-95085526(+)', 'TCACAGTGGCTAAGTTCTGCA')]
        obs = read_fasta(self.fp_fasta_mirna)
        self.assertEqual(exp, list(obs))

    def test_read_CT_file(self):
        with self.assertRaises(ValueError) as e:
            obs = list(read_CT_file(getFP('data/broken.ct')))
        with self.assertRaises(ValueError) as e:
            obs = list(read_CT_file(getFP('data/broken_missing.ct')))

        exp = [
            ('A stem-loop structure (with a bulge)',
             'GGGCAAUCCUCUUCGGGCCC',
             {1: 20, 2: 19, 3: 18, 4: 17, 8: 16, 9: 15, 15: 9, 16: 8, 17: 4, 18: 3, 19: 2, 20: 1}),
            ('A pseudo-knot structure',
             'GAUGGCACUCCCAUCAAUUGGAGC',
             {1: 15, 2: 14, 3: 13, 4: 12, 5: 11, 8: 23, 9: 22, 10: 21, 11: 5, 12: 4, 13: 3, 14: 2, 15: 1, 21: 10, 22: 9, 23: 8})]
        obs = read_CT_file(self.fp_multiCT)

        obs = list(read_CT_file(getFP('data/SARS-CoV-2_Full_Length_Secondary_Structure_Map.ct')))[0]
        self.assertEqual('finalStructure wPKs', obs[0])
        self.assertEqual(29903, len(obs[1]))
        self.assertTrue(obs[1].endswith('ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'))
        self.assertTrue(len(obs[1]) > len(obs[2]))

    def test_disentangle_knots(self):
        obs = disentangle_knots(next(read_CT_file(getFP('data/earlyKnot.ct')))[2], verbose=None)
        self.assertEqual(len(obs['nested']), 6*2)
        self.assertEqual(len(obs['knotted']), 3*2)

        obs = disentangle_knots(next(read_CT_file(getFP('data/lateKnot.ct')))[2], verbose=None)
        self.assertEqual(len(obs['nested']), 6*2)
        self.assertEqual(len(obs['knotted']), 3*2)

        obs = disentangle_knots(next(read_CT_file(self.fp_sars))[2], verbose=None)
        self.assertEqual(len(obs['nested']), 8681*2)
        self.assertEqual(len(obs['knotted']), 7*2)

if __name__ == '__main__':
    main()
