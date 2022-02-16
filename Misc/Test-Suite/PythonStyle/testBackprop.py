from unittest import TestCase, main

from test_utils import testBackprop

class BackpropagationTest(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_deadends(self):
        # issue with NTs visited intermediately but not used in the end,
        # e.g. ml_comps1 with cadd(x,y) where x might be a valid parse, but not y
        inpseq = 'CCCaaaCaaaGG'
        self.assertFalse(testBackprop(inpseq))

    def test_firstexample(self):
        # first input sequence that worked for pfunc
        inpseq = 'GGaaaCaaaCCC'
        self.assertFalse(testBackprop(inpseq))

    def test_deadendsInMultiloops(self):
        # small RNA with multiloops
        # I here saw that any occurence in the param array should lead to termination
        inpseq = 'CCCCaaaGGCCaaaGGGG'
        self.assertFalse(testBackprop(inpseq))

    def test_largerExample(self):
        # max size example, otherwise enum explodes
        inpseq = 'ACCCUACUGUGCUAACCGAACCAGA'
        self.assertFalse(testBackprop(inpseq))

if __name__ == '__main__':
    main()
