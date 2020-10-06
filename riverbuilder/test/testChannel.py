import unittest
import numpy as np
from math import pi, inf, sqrt
import matplotlib.pyplot as plt

from ..core.cChannel import Channel

class CoreChannelTest(unittest.TestCase):

    def test_cutArea(self):
        '''
        '''
        avail_pts = [[{1,2,3,4,5}],
                    [set(list(range(10)))],
                    [set(list(range(10)))],
                    [set(list(range(10)))],
                    [set(list(range(10)))],
                    [set(list(range(10)))],
                    [{1,2,3,4,5}]
                    ]
        size = 3
        check = set(list(range(7)))
        xmin = 0
        xmax = 7
        c = Channel()
        area, check = c.cutArea(avail_pts, size, check, xmin, xmax)
#        print('area', area)
#        print('check', check)
#        print('avail_pts', avail_pts)

        size = 15
        area, check = c.cutArea(avail_pts, size, check, xmin, xmax)
        self.assertEqual(area, [])

    def test_createBoulder(self):
        c = Channel()
        area = [(0, 10), (0, 10)]
        b = c.createBoulder(area)
        for i in b:
            print(i)

        area = [(-4, 5), (-4, 5)]
        b = c.createBoulder(area)
        for i in b:
            print(i)
        
if __name__ == '__main__':
    unittest.main()
