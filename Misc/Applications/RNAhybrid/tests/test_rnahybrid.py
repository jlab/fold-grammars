import sys
from os.path import dirname, join
sys.path.append(dirname(__file__) + '/../')

from unittest import TestCase, main
from input import read_CT_file, disentangle_knots

def getFP(filepath):
    return join(dirname(__file__), filepath)

class TestRNAhybrid(TestCase):
    def setUp(self):
        self.CTsars = disentangle_knots(next(read_CT_file(getFP('data/SARS-CoV-2_Full_Length_Secondary_Structure_Map.ct')))[2], verbose=None)['nested']
        # Pairs are stored as opening: closing base pair positions, i.e. opening < closing.
        # To ease access, we enrich this dict by also adding in closing: opening information.
        self.CTsars.update({c: o for o, c in self.CTsars.items()})

    def tearDown(self):
        pass

    def test_getOriginalPairs(self):
        exp = {839: 851, 840: 850, 841: 849, 842: 848, 854: 943, 855: 942, 856: 941, 857: 940, 862: 906,  # assert that: i < j
               848: 842, 850: 840, 852: 838, 853: 837  # duplicate pairs where j < i
              }

        novel_target_pairing_partners = [int(x) for x in '839,840,841,842,843,845,846,847,848,850,852,853,854,855,856,857,858,859,862'.split(',')]
        obs = {o: c for o, c in {x: self.CTsars.get(x, -1) for x in novel_target_pairing_partners}.items() if c != -1}
        self.assertEqual(exp, obs)


if __name__ == '__main__':
    main()
