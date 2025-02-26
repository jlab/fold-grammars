import sys
from os.path import dirname
sys.path.append(dirname(__file__) + '/../')

from unittest import TestCase, main

from parse import *
from execute import read_CT_file, disentangle_knots

class TestExecute(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read_CT_file(self):
        obs = read_CT_file('tests/data/SARS-CoV-2_Full_Length_Secondary_Structure_Map.ct')
        self.assertEqual('wPKs', obs['name'])
        self.assertEqual(29903, len(obs['sequence']))
        self.assertTrue(obs['sequence'].endswith('ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'))

        with self.assertRaises(ValueError) as e:
            obs = read_CT_file('tests/data/broken.ct')

        with self.assertRaises(ValueError) as e:
            obs = read_CT_file('tests/data/broken_missing.ct')

    def test_disentangle_knots(self):
        obs = disentangle_knots(read_CT_file('tests/data/earlyKnot.ct')['structure'], verbose=None)
        self.assertEqual(len(obs['nested']), 6)
        self.assertEqual(len(obs['knotted']), 3)

        obs = disentangle_knots(read_CT_file('tests/data/lateKnot.ct')['structure'], verbose=None)
        self.assertEqual(len(obs['nested']), 6)
        self.assertEqual(len(obs['knotted']), 3)

if __name__ == '__main__':
    main()
