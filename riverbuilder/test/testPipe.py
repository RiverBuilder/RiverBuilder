import unittest
import numpy as np
from math import pi, inf, sqrt
import matplotlib.pyplot as plt
from datetime import datetime

from ..core.cPipe import Pipe

class CoreChannelTest(unittest.TestCase):

    def test_setCenterline(self):
        ''' Case 1: no function involved
            create straight horizontal centerline
        '''
        x_len = 100
        x_slope = 0
        pipe = Pipe(x_len, x_slope)
        pipe.setCenterline()

        xCenterExpected = np.array(list(range(100)))
        yCenterExpected = np.array([0]*100)
        self.assertTrue(np.array_equal(xCenterExpected, pipe.x_v))
        self.assertTrue(np.array_equal(yCenterExpected, pipe.y_center))

        ''' Case 2: fun: f(x) = 2*x was applied to centerline
            create straight horizontal centerline
        '''
        fun = lambda x : (x, 1.5*x)
        pipe.setCenterline(fun)

        xCenterExpected = np.array(list(range(100)))
        yCenterExpected = np.array(1.5*xCenterExpected)
        self.assertTrue(np.array_equal(xCenterExpected, pipe.x_v))
        self.assertTrue(np.array_equal(yCenterExpected, pipe.y_center))

    def test_setCenterlineManual(self):
        ''' Case 1: given x_v, y_v with equal length.
            Pipe's x_v, y_center should be replaced.
        '''
        x_len = 100
        x_slope = 0
        pipe = Pipe(x_len, x_slope)
        pipe.setCenterline()

        xCenterExpected = np.array(list(range(50)))
        yCenterExpected = np.array(list(range(50)))
        pipe.setCenterlineManual(yCenterExpected, xCenterExpected)
        self.assertTrue(np.array_equal(xCenterExpected, pipe.x_v))
        self.assertTrue(np.array_equal(yCenterExpected, pipe.y_center))

        ''' Case 2: x_v is skipped, y_v has the same length as pipe's original y_center
            Only pipe's y_center should be replaced.
        '''
        pipe = Pipe(x_len, x_slope)
        pipe.setCenterline()

        xCenterExpected = np.array(list(range(100)))
        yCenterExpected = np.array(list(range(100)))
        pipe.setCenterlineManual(yCenterExpected)
        self.assertTrue(np.array_equal(xCenterExpected, pipe.x_v))
        self.assertTrue(np.array_equal(yCenterExpected, pipe.y_center))

        ''' Case 3: x_v is skipped, y_v is longer than pipe's original y_center
            Pipe's y_center should be replaced with respect to the length of x_v
        '''
        pipe = Pipe(x_len, x_slope)
        pipe.setCenterline()

        xCenterExpected = np.array(list(range(100)))
        yCenterExpected = np.array(list(range(100)))
        y_v = np.array(list(range(150)))
        pipe.setCenterlineManual(y_v)
        self.assertTrue(np.array_equal(xCenterExpected, pipe.x_v))
        self.assertTrue(np.array_equal(yCenterExpected, pipe.y_center))

        ''' Case 4: x_v is skipped, y_v is shorter than pipe's original y_center
            Pipe's y_center should be replaced with respect to the length of x_v
        '''
        pipe = Pipe(x_len, x_slope)
        pipe.setCenterline()

        xCenterExpected = np.array(list(range(100)))
        yCenterExpected = np.array(list(range(50))+[49]*50)
        y_v = np.array(list(range(50)))
        pipe.setCenterlineManual(y_v)
        self.assertTrue(np.array_equal(yCenterExpected, pipe.y_center))

    def test_smoothCenterline(self):
        ''' Case 1: smooth degree = 0; no smoothen at all.
            Should get original centerline.
        '''
        x_len = 10
        x_slope = 0
        pipe = Pipe(x_len, x_slope)
        y_v = np.array(list(range(10)))
        pipe.setCenterlineManual(y_v)

        yExpected = y_v
        pipe.smoothCenterline(0)
        self.assertTrue(np.array_equal(yExpected, pipe.y_center))

        ''' Case 2: smooth degree = 1; take average with every 3 points.
        '''
        pipe = Pipe(x_len, x_slope)
        y_v = np.array([0, 6, 3, 0, -3, -6, -3, 0, 3, 6])
        pipe.setCenterlineManual(y_v)

        yExpected = np.array([3, 3, 3, 0, -3, -4, -3, 0, 3, 4.5])
        pipe.smoothCenterline(1)
        self.assertTrue(np.array_equal(yExpected, pipe.y_center))

    def test_addTopOffset(self):
        '''
        A simple case
        '''
        newLevel_x = np.array(list(range(10)))
        prevLevel_x = np.array(list(range(10)))
        newLevel_y = np.array([1]*10)
        prevLevel_y = np.array([0]*10)
        prevLevel_z = np.array(list(range(10)))
        offset = 1
        pipe = Pipe()

        expected = prevLevel_z + 1
        got = pipe.addTopOffset(newLevel_x, newLevel_y, prevLevel_x, prevLevel_y, prevLevel_z, offset)
        self.assertTrue(np.array_equal(expected, got))

        '''
        Test the speed
        '''
        newLevel_x = np.array(list(range(10000)))
        prevLevel_x = np.array(list(range(10000)))
        newLevel_y = np.array([1]*10000)
        prevLevel_y = np.array([0]*10000)
        prevLevel_z = np.array(list(range(10000)))

        expected = prevLevel_z + 1
        time1 = datetime.now()
        got = pipe.addTopOffset(newLevel_x, newLevel_y, prevLevel_x, prevLevel_y, prevLevel_z, offset)
        self.assertTrue(np.array_equal(expected, got))
        time2 = datetime.now()
        print('delta time is:', time2-time1)

        '''
        A more complicated cases.
        '''
        newLevel_x = np.array([1, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        prevLevel_x = np.array(list(range(10)))
        newLevel_y = np.array([1]*10)
        prevLevel_y = np.array([0]*10)
        prevLevel_z = np.array(list(range(10)))
        offset = 1
        pipe = Pipe()

        expected = np.copy(newLevel_x) + 1
        got = pipe.addTopOffset(newLevel_x, newLevel_y, prevLevel_x, prevLevel_y, prevLevel_z, offset)
        self.assertTrue(np.array_equal(expected, got))

        '''
        When newLevel has fewer points then previous level.
        '''
        newLevel_x = np.array([0])
        newLevel_y = np.array([1])
        expected = np.array([1])
        got = pipe.addTopOffset(newLevel_x, newLevel_y, prevLevel_x, prevLevel_y, prevLevel_z, offset)
        self.assertTrue(np.array_equal(expected, got))

        '''
        When newLevel has more points then previous level.
        '''
        newLevel_x = np.array(list(range(10)))
        prevLevel_x = np.array([2, 7])
        newLevel_y = np.array([1]*10)
        prevLevel_y = np.array([0, 0])
        prevLevel_z = np.array([1, 3])
        offset = 1
        pipe = Pipe()

        expected = np.array([2]*5+[4]*5)
        got = pipe.addTopOffset(newLevel_x, newLevel_y, prevLevel_x, prevLevel_y, prevLevel_z, offset)
        self.assertTrue(np.array_equal(expected, got))

    def test_getDynamicPipeSlope(self):
        '''
        A simple case where the pipe is horizontal, 
        therefore, dynamic pipe slope should be the repeat of x slope
        '''
        pipe = Pipe(x_len=10)
        expected = np.array([0.01] * 10)
        got = pipe.getDynamicPipeSlope()
        self.assertTrue(np.array_equal(expected, got))

        '''
        case where the first half of the pipe is horizontal,
        while the second half of the pipe has sinuosity of 2.
        '''
        pipe.s_center = np.array(list(range(5)) + [6, 8, 10, 12, 14])
        expected = np.array([0.01]*5+[0.005]*5)
        got = pipe.getDynamicPipeSlope()
        self.assertTrue(np.array_equal(expected, got))


if __name__ == '__main__':
    unittest.main()
