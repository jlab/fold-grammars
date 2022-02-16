from unittest import TestCase, main

from test_utils import testBackprop

class BackpropagationTest(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_backprop(self):
        # issue with NTs visited intermediately but not used in the end,
        # e.g. ml_comps1 with cadd(x,y) where x might be a valid parse, but not y
        inpseq = 'CCCaaaCaaaGG'
        self.assertFalse(testBackprop(inpseq))

        # first input sequence that worked for pfunc
        inpseq = 'GGaaaCaaaCCC'
        self.assertFalse(testBackprop(inpseq))

        # small RNA with multiloops
        inpseq = 'CCCCaaaGGCCaaaGGGG'
        self.assertFalse(testBackprop(inpseq))

if __name__ == '__main__':
    main()
