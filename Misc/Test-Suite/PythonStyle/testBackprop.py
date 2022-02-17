from unittest import TestCase, main

from test_utils import testBackprop, testBasepair
import nodangle as nd
import pandas as pd

class BackpropagationTest(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_deadends(self):
        # issue with NTs visited intermediately but not used in the end,
        # e.g. ml_comps1 with cadd(x,y) where x might be a valid parse, but not y
        inpseq = 'CCCaaaCaaaGG'
        self.assertFalse(testBackprop(inpseq, verbose=None))

    def test_firstexample(self):
        # first input sequence that worked for pfunc
        inpseq = 'GGaaaCaaaCCC'
        self.assertFalse(testBackprop(inpseq, verbose=None))

    def test_deadendsInMultiloops(self):
        # small RNA with multiloops
        # I here saw that any occurence in the param array should lead to termination
        inpseq = 'CCCCaaaGGCCaaaGGGG'
        self.assertFalse(testBackprop(inpseq, verbose=None))

    def test_largerExample(self):
        # max size example, otherwise enum explodes
        inpseq = 'ACCCUACUGUGCUAACCGAACCAGA'
        self.assertFalse(testBackprop(inpseq, verbose=None))

    def test_switchAlgebra(self):
        # issue with NTs visited intermediately but not used in the end,
        # e.g. ml_comps1 with cadd(x,y) where x might be a valid parse, but not y
        inpseq = 'CCCaaaCaaaGG'
        self.assertFalse(testBackprop(inpseq, algebra='count', verbose=None))

class BasepairTest(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_small_isolated(self):
        inpseq = 'ACCCUACUGUGCUAACCGAACCAGA'
        self.assertFalse(testBasepair(inpseq, verbose=None))

    def test_A_isolated(self):
        inpseq = 'AAGGGCGUCGUCGCCCCGAGUCGUAGCAGUUGACUACUGUUAUGU'
        self.assertFalse(testBasepair(inpseq, verbose=None))

    def test_B_isolated(self):
        inpseq = 'gGGCCGGGCGCGGUGGCGCGCGCCUGUAGUCCCAGCUACUCGGGAGGCUC'
        self.assertFalse(testBasepair(inpseq, verbose=None))

    def test_C_isolated(self):
        inpseq = 'ACCCUACUGUGCUAACCGAACCAGAUAACGGUACAGUAGGGGUAAAUUCUCCGCAUUCGGUGCGGAAAA'
        self.assertFalse(testBasepair(inpseq, verbose=None))

class MFETest(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_mfe(self):
        inpseq = 'CCCCaaaGGCCaaaGGGG'
        nd.init(inpseq, algebra='mfe', printstack=False, printBTstack=False)
        fwd = nd.nt_struct(0)

        # this means, that backprop is also good for non multiplicative algebras :-)
        exp = -330  # manually determined
        obs = nd.backtrace(3,17,'leftB')
        self.assertEqual(exp, obs)

class FwdTabulationTest(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_fwd(self):
        res = []
        for design in [
            ["hairpin","leftB","multiloop","rightB","stack","dangle","iloop","ml_comps","ml_comps1","strong","struct","weak"],
            ["dangle","iloop","ml_comps","ml_comps1","strong","struct","weak"],
            ["dangle","iloop","ml_comps","ml_comps1","strong","struct"],
            ["iloop","ml_comps","ml_comps1","strong","struct"],
            ["ml_comps","ml_comps1","struct"],
            ["ml_comps","struct"],
            ["struct"],
            [],
            ]:
            nd.init('ACCCUACUGUGCUAACCGAACCAGA', algebra='pfunc', printstack=False, printBTstack=False, tabulateNTs=design)
            res.append({
                'result': nd.nt_struct(0),
                'steps': nd.COMPUTATIONALSTEPS,
                'storage': nd.STORAGE,
                'design': design})
        res = pd.DataFrame(res)

        self.assertEqual(list(res['result'].unique()), [0.0008308523778764974])
        self.assertEqual(list(res['steps'].values), [5018, 5018, 5540, 5961, 126793, 180479, 278748, 745656])
        self.assertEqual(list(res['storage'].values), [1922, 1054, 823, 592, 130, 78, 26, 0])


if __name__ == '__main__':
    main()
