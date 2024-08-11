import sys
from os.path import dirname
sys.path.append(dirname(__file__) + '/../')

from unittest import TestCase, main

from parse import *

class TestParse(TestCase):
    def setUp(self):
        self.t1 = Product(TypeInt('a'))
        self.t2 = Product(TypeInt('a'), TypeInt('b'))
        self.t3 = Product(TypeInt('a'), Product(TypeInt('b'), TypeInt('c')))
        self.t4 = Product(Product(TypeInt('a'), TypeInt('b')), TypeInt('c'))

    def tearDown(self):
        pass

    def test_product_to_pattern(self):
        exp = r'(-?\d+)'
        self.assertEqual(self.t1.getRegex(), exp)

        exp = r'\( (-?\d+) , (-?\d+) \)'
        self.assertEqual(self.t2.getRegex(), exp)

        exp = r'\( (-?\d+) , \( (-?\d+) , (-?\d+) \) \)'
        self.assertEqual(self.t3.getRegex(), exp)

        exp = r'\( \( (-?\d+) , (-?\d+) \) , (-?\d+) \)'
        self.assertEqual(self.t4.getRegex(), exp)

    def test_types(self):
        m = TypeMFE()
        self.assertEqual(m.getRegex(), r'(-?\d+\.?\d*)')

        with self.assertRaises(ValueError) as e:
            p = Product(self.t2)
        self.assertTrue('Cannot make a Product of a single Product!' in str(e.exception))

    def test_parse_lines(self):
        obs = self.t2.parse_lines(["( 23 , -10 )", "( -44 , 10 )", "( -23 , -35347547568 )"])
        self.assertEqual(obs, [
            {'a': 23, 'b': -10},
            {'a': -44, 'b': 10},
            {'a': -23, 'b': -35347547568}])

        p = Product(TypeHybrid())
        exp = [{'target_position': 25682,
                'target_unpaired': 'U     C       GA      U U',
                'target_stacked': ' UGCUG UGGCCUU  AGCCCC U ',
                'pairs': ' ||||| |:|||:|  |||||| : ',
                'mirna_stacked': ' ACGAC AUCGGGA--UCGGGG G ',
                'mirna_unpaired': '      A               C U'}]
        obs = p.parse_lines(["(25682, target 5' U     C       GA      U U 3',            UGCUG UGGCCUU  AGCCCC U ,            ||||| |:|||:|  |||||| : ,            ACGAC AUCGGGA--UCGGGG G , miRNA 3'        A               C U 5')"])
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
