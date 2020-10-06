import unittest
import numpy as np
from math import pi, inf, sqrt

from ..generators import functionTransformer

class FunctionTransformerTests(unittest.TestCase):

    def test_transformTri(self):
        # reah_len is the same as tread_len
        # ExpectOut: "(1,1,0,MASK2)"
        function = [1, 1, 0, 'MASK1']
        maskFunc = "MASK2"
        reach_len = 100
        tread_len = 100
        startPoint = (0, 0)
        yExpect = "(1, 1.0, 0, MASK2)"
        yOut = functionTransformer.transformTri(function, maskFunc, reach_len, tread_len, startPoint)
        self.assertEqual(yExpect, yOut)

        # reah_len is the half as tread_len
        # ExpectOut: "(1,1,0,MASK2)"
        function = [1, 1, 0, 'MASK1']
        maskFunc = "MASK2"
        reach_len = 100
        tread_len = 50
        startPoint = (0, 0)
        yExpect = "(1, 2.0, 0, MASK2)"
        yOut = functionTransformer.transformTri(function, maskFunc, reach_len, tread_len, startPoint)
        self.assertEqual(yExpect, yOut)

#        # startPoint is at middle
#        # ExpectOut: "(1,1,0,MASK2)"
#        function = [1, 1, 0, 'MASK1']
#        maskFunc = "MASK2"
#        reach_len = 100
#        tread_len = 50
#        startPoint = (50, 0)
#        yExpect = "(1, 2.0, 0, MASK2)"
#        yOut = functionTransformer.transformTri(function, maskFunc, reach_len, tread_len, startPoint)
#        print(yExpect)
#        print(yOut)
#        self.assertTrue(np.array_equal(yExpect, yOut))
#
#        # startPoint is at quarter
#        # ExpectOut: "(1,1,0,MASK2)"
#        function = [1, 1, 0, 'MASK1']
#        maskFunc = "MASK2"
#        reach_len = 100
#        tread_len = 50
#        startPoint = (0, 0)
#        yExpect = "(1, 2.0, pi, MASK2)"
#        yOut = functionTransformer.transformTri(function, maskFunc, reach_len, tread_len, startPoint)
#        #print(yExpect)
#        #print(yOut)
#        self.assertTrue(np.array_equal(yExpect, yOut))

    def test_transform(self):
        # Invalid case
        # ExpectOut: ''
        funcName = 'XX1'
        function = [1, 1, 0, 'MASK1']
        maskFunc = "MASK2"
        reach_len = 100
        tread_len = 100
        startPoint = (0, 0)
        yExpect = ''
        yOut = functionTransformer.transform(funcName, function, maskFunc, reach_len, tread_len, startPoint)
        self.assertEqual(yExpect, yOut)

        # A very simple case
        # ExpectOut: "(1,1,0,MASK2)"
        funcName = 'SIN13'
        function = [1, 1, 0, 'MASK1']
        maskFunc = "MASK2"
        reach_len = 100
        tread_len = 100
        startPoint = (0, 0)
        yExpect = "SIN=(1, 1.0, 0, MASK2)"
        yOut = functionTransformer.transform(funcName, function, maskFunc, reach_len, tread_len, startPoint)
        self.assertEqual(yExpect, yOut)


if __name__ == '__main__':
    unittest.main()
