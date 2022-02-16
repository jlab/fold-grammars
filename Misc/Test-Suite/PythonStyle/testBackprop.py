from unittest import TestCase, main

from test_utils import testBackprop

class BackpropagationTest(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_backprop(self):
        self.assertFalse(testBackprop('CCCaaaCaaaGG'))


if __name__ == '__main__':
    main()
