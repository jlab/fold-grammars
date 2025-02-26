import sys
from os.path import dirname, join
sys.path.append(dirname(__file__) + '/../')

from unittest import TestCase, main

from parse_gapc import *
from input import read_CT_file, disentangle_knots

def getFP(filepath):
    return join(dirname(__file__), filepath)

class TestExecute(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass



if __name__ == '__main__':
    main()
