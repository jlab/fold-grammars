import sys
from os.path import dirname, join
sys.path.append(dirname(__file__) + '/../')

from unittest import TestCase, main
from execute import read_CT_file, disentangle_knots
from RNAhybrid import extend_pairs_to_valid_substructure, get_sub_dotBracket_structure

def getFP(filepath):
    return join(dirname(__file__), filepath)

class TestRNAhybrid(TestCase):
    def setUp(self):
        self.CTsars = disentangle_knots(read_CT_file(getFP('data/SARS-CoV-2_Full_Length_Secondary_Structure_Map.ct'))['structure'], verbose=None)['nested']
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

    def test_extend_pairs_to_valid_substructure(self):
        exp_db = '((((((.....))))))((((..((((((((..(((..((((((((....)))))..)))))).)))))))).((((........))))..............))))' # 31 pairs

        exp_pairs = {
            839: 851, 840: 850, 841: 849, 842: 848, 848: 842, 850: 840,
            852: 838, 853: 837, 854: 943, 855: 942, 856: 941, 857: 940,
            862: 906, 837: 853, 838: 852, 851: 839, 849: 841, 943: 854,
            942: 855, 941: 856, 940: 857, 860: 908, 908: 860, 861: 907,
            907: 861, 906: 862, 863: 905, 905: 863, 864: 904, 904: 864,
            865: 903, 903: 865, 866: 902, 902: 866, 867: 901, 901: 867,
            870: 899, 899: 870, 871: 898, 898: 871, 872: 897, 897: 872,
            875: 896, 896: 875, 876: 895, 895: 876, 877: 894, 894: 877,
            878: 891, 891: 878, 879: 890, 890: 879, 880: 889, 889: 880,
            881: 888, 888: 881, 882: 887, 887: 882, 910: 925, 925: 910,
            911: 924, 924: 911, 912: 923, 923: 912, 913: 922, 922: 913}

        novel_target_pairing_partners = [int(x) for x in '839,840,841,842,843,845,846,847,848,850,852,853,854,855,856,857,858,859,862'.split(',')]
        original_target_pairs = {o: c for o, c in {x: self.CTsars.get(x, -1) for x in novel_target_pairing_partners}.items() if c != -1}

        obs_ext_pairs = extend_pairs_to_valid_substructure(self.CTsars, original_target_pairs)
        self.assertEqual(exp_pairs, obs_ext_pairs)

        h = sorted(list(obs_ext_pairs.keys()) + list(obs_ext_pairs.values()))
        obs = get_sub_dotBracket_structure(self.CTsars, h[0], h[-1])
        self.assertEqual(exp_db, obs)

        pair_subset = {28395: 28379, 28396: 28377, 28397: 28376, 28398: 28375, 28399: 28374, 28410: 28404, 28411: 28403}
        exp = {28395: 28379, 28396: 28377, 28397: 28376, 28398: 28375, 28399: 28374, 28410: 28404, 28411: 28403, 28374: 28399, 28375: 28398, 28376: 28397, 28377: 28396, 28379: 28395, 28380: 28394, 28394: 28380, 28381: 28393, 28393: 28381, 28382: 28392, 28392: 28382, 28400: 28373, 28373: 28400, 28401: 28372, 28372: 28401, 28403: 28411, 28404: 28410}
        obs = extend_pairs_to_valid_substructure(self.CTsars, pair_subset)
        self.assertEqual(exp, obs)

if __name__ == '__main__':
    main()
