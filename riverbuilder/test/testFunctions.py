import unittest
import numpy as np
from math import pi, inf, sqrt
import matplotlib.pyplot as plt

from ..core import functions

class CoreFunctionsTest(unittest.TestCase):

    def test_maskCalculate(self):
        '''
        A simple test where mask is applied to all points.
        '''
        x_v = np.array(list(range(20)))
        y_v = np.array(list(range(20)))
        mask = [[0, 20, True]]
        x_expected = np.array(list(range(20)))
        y_expected = np.array(list(range(20)))
        x_got, y_got = functions.maskCalculate(y_v, x_v, mask)
        self.assertTrue(np.array_equal(x_expected, x_got))
        self.assertTrue(np.array_equal(y_expected, y_got))

        '''
        A simple test where mask is applied to no points.
        '''
        mask = [[0, 10, False], [10, 20, False]]
        x_expected = np.array(list(range(20)))
        y_expected = np.array([0]*20)
        x_got, y_got = functions.maskCalculate(y_v, x_v, mask)
        self.assertTrue(np.array_equal(x_expected, x_got))
        self.assertTrue(np.array_equal(y_expected, y_got))

        '''
        A test where mask is applied to first half.
        '''
        mask = [[0, 10, True], [10, 20, False]]
        x_expected = np.array(list(range(20)))
        y_expected = np.array(list(range(10))+[0]*10)
        x_got, y_got = functions.maskCalculate(y_v, x_v, mask)
        self.assertTrue(np.array_equal(x_expected, x_got))
        self.assertTrue(np.array_equal(y_expected, y_got))

    def test_deleteCycle(self):
        line = [(0, 3), (1, 2), (2, 1), (3, 0), (4, -1), (5, -2), (4.5, 2.5), (2.8, -0.2), (1.5, -1)]
        expected = [(0, 3), (1, 2), (2, 1), (2.8, -0.2), (1.5, -1)]
        out = functions.deleteCycles(line)
#        print(out)
        self.assertEqual(expected, out)

        line = [(0, 0), (1, 0), (2, 0),(1.6, 1), (1.5, 1), (1.4, 1), (1.4, -1)]
        expected = [(0, 0), (1, 0), (1.4, -1)]
        out = functions.deleteCycles(line)
#        print(out)
        self.assertEqual(expected, out)
if __name__ == '__main__':
    unittest.main()
