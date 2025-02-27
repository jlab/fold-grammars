import sys
from os.path import dirname, join
sys.path.append(dirname(__file__) + '/../')

from unittest import TestCase, main

from input import read_fasta, read_CT_file, disentangle_knots, nested_pairs_to_dotBracket, get_minimal_valid_substructure

def getFP(filepath):
    return join(dirname(__file__), filepath)

class TestExecute(TestCase):
    def setUp(self):
        self.fp_fasta_mirna = getFP('data/mirnas.fasta')
        self.fp_multiCT = getFP('data/multi.ct')
        self.fp_sars = getFP('data/SARS-CoV-2_Full_Length_Secondary_Structure_Map.ct')
        self.pairs1 = {14292: 14269, 14293: 14268, 14307: 14437, 14308: 14436, 14309: 14323}
        self.pairs2 = {27149: 27127, 27150: 27126, 27164: 27330, 27165: 27329}
        self.pairs3 = {28395: 28379, 28396: 28377, 28397: 28376, 28398: 28375, 28399: 28374, 28410: 28404, 28411: 28403}
        self.pairs4 = {4060: 4053, 4061: 4052, 4062: 4051, 4063: 4050, 4067: 4048, 4068: 4047, 4069: 4046, 4071: 4144, 4072: 4143, 4073: 4142, 4074: 4141, 4075: 4140, 4076: 4139}

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

    def test_nested_pairs_to_dotBracket(self):
        exp = '((......................)).............(((.............)................................................................................................................))'
        obs = nested_pairs_to_dotBracket(self.pairs1)
        self.assertEqual(exp, obs)

        exp = '((.....................)).............((...................................................................................................................................................................))'
        obs = nested_pairs_to_dotBracket(self.pairs2)
        self.assertEqual(exp, obs)

        exp = '((((.(...............)))))...((.....))'
        obs = nested_pairs_to_dotBracket(self.pairs3)
        self.assertEqual(exp, obs)

        exp = '(((.((((......))))...))).((((((..............................................................))))))'
        obs = nested_pairs_to_dotBracket(self.pairs4)
        self.assertEqual(exp, obs)

        exp = '....(....)..............((...(((....)))..))...............'
        obs = nested_pairs_to_dotBracket(
            {24728: 24741, 24729: 24740, 24730: 24739, 24731: 24738, 24732: 24737, 24742: 24785, 24743: 24784, 24744: 24783, 24745: 24781, 24746: 24780, 24747: 24779, 24748: 24776, 24749: 24775, 24750: 24772, 24751: 24771, 24752: 24770, 24753: 24769, 24757: 24766, 24758: 24765, 24759: 24764},
            [24738, 24739, 24740, 24741, 24742, 24743, 24744, 24745, 24746, 24747, 24748, 24749, 24750, 24751])
        self.assertEqual(exp, obs)


    def test_get_minimal_valid_substructure(self):
        fullSARS = list(read_CT_file(getFP('data/SARS-CoV-2_Full_Length_Secondary_Structure_Map.ct')))[0]

        obs = get_minimal_valid_substructure(fullSARS[2], self.pairs1)
        self.assertEqual(nested_pairs_to_dotBracket(obs), '((((.((((((.......))))))..))))......(((((((((...........))..((((((.((....)).))))))..(...(((((.....)))))...).(((......)))...(((((.((((((.(((..............))).)))).))))))).)))))))')

        obs = get_minimal_valid_substructure(fullSARS[2], self.pairs2)
        self.assertEqual(nested_pairs_to_dotBracket(obs), '(((((.(((.(((.....................))).))).)))))((((((...(((((((((((.((((...((((..(((((((((...))))).))))))))....).)))...))))))).......)))).(((((((((.......(((.((((....)))).)))..............)))))))))...............))))))')

        obs = get_minimal_valid_substructure(fullSARS[2], self.pairs3)
        self.assertEqual(nested_pairs_to_dotBracket(obs), '((((((.((((.........)))))))))).((.....))')

        obs = get_minimal_valid_substructure(fullSARS[2], self.pairs4)
        self.assertEqual(nested_pairs_to_dotBracket(obs), '(((((((((......))))..)))))((((((.....(((...((((((((((.((......................)).)))))))))))))))))))')

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
        original_target_pairs = {o: c for o, c in {x: fullSARS[2].get(x, -1) for x in novel_target_pairing_partners}.items() if c != -1}

        obs_ext_pairs = get_minimal_valid_substructure(fullSARS[2], original_target_pairs)
        self.assertEqual(exp_pairs, obs_ext_pairs)

        obs = nested_pairs_to_dotBracket(obs_ext_pairs)
        self.assertEqual(exp_db, obs)

        pair_subset = {28395: 28379, 28396: 28377, 28397: 28376, 28398: 28375, 28399: 28374, 28410: 28404, 28411: 28403}
        exp = {28395: 28379, 28396: 28377, 28397: 28376, 28398: 28375, 28399: 28374, 28410: 28404, 28411: 28403, 28374: 28399, 28375: 28398, 28376: 28397, 28377: 28396, 28379: 28395, 28380: 28394, 28394: 28380, 28381: 28393, 28393: 28381, 28382: 28392, 28392: 28382, 28400: 28373, 28373: 28400, 28401: 28372, 28372: 28401, 28403: 28411, 28404: 28410}
        obs = get_minimal_valid_substructure(fullSARS[2], pair_subset)
        self.assertEqual(exp, obs)

if __name__ == '__main__':
    main()
