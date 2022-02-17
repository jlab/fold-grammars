from unittest import TestCase, main

from test_utils import testBackprop, testBasepair
import nodangle as nd

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

if __name__ == '__main__':
    main()
