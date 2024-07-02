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

    # def test_parse_lines(self):
    #     obs = Product(self.t2).parse_lines(["( 23 , -10 )"])
    #
    #     obs = Product(self.t1).parse_lines(["-23"])
        #print(obs)
        #obs = parse_gapc(["23", "-23", "-35347547568"], self.t1)
        # self.assertEqual(obs, [[23], [-23], [-35347547568]])

if __name__ == '__main__':
    main()
