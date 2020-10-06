import unittest
import numpy as np
from math import pi, inf, sqrt

from ..generators import functionGenerator

class FunctionGeneratorTest(unittest.TestCase):

    def test_intToStrComplete(self):
        '''A very simple case
        outExpected - '1'
        '''
        x = 1
        xtot = 1
        outExpected = '1'
        outGot = functionGenerator.intToStrComplete(x, xtot)
        self.assertEqual(outExpected, outGot)

        '''A case need completion
        outExpected - '001'
        '''
        x = 1
        xtot = 123
        outExpected = '001'
        outGot = functionGenerator.intToStrComplete(x, xtot)
        self.assertEqual(outExpected, outGot)


    def test_generateMask(self):
        '''A very simple case
        outExpected - 'MASK1=(0, on, 50, off)'
        '''
        i = '1'
        start = '0'
        end = '50'
        outExpected = 'MASK1=(0, on, 50, off)'
        outGot = functionGenerator.generateMask(i, start, end)
        self.assertEqual(outExpected, outGot)

        '''Another case
        outExpected - 'MASK10203=(20, on, 50, off)'
        '''
        i = '10203'
        start = '20'
        end = '50'
        outExpected = 'MASK10203=(20, on, 50, off)'
        outGot = functionGenerator.generateMask(i, start, end)
        self.assertEqual(outExpected, outGot)
        
        
        #def test_parseFile(self):
        #'''
        #Input file is inputSamples/stepPoolThalwegBasic.txt
        #'''
        #path = 'inputSamples/stepPoolThalwegBasic.txt'
        #reachPara, stepsList = functionGenerator.parseFile(path)
        ##print(reachPara)
        ##for step in stepsList:
        ##    print(step)
        #check = True
        #check = check and (reachPara['Length'] == 100)
        #check = check and (reachPara['Slope'] == 0.0)
        #check = check and (reachPara['MASK0'] == 'ALL')
        #check = check and (reachPara['SIN1'] == '(1, 1, 0, MASK0)')
        #check = check and (reachPara['COS2'] == '(1, 1, 0, MASK0)')
        #check = check and (stepsList[0]['step function'] == "flat")
        #check = check and (stepsList[1]['step length'] == 10)
        #check = check and (stepsList[1]['step height'] == 2)
        #check = check and (stepsList[1]['connector length'] == 1)
        #check = check and (stepsList[1]['step function'] == "SIN1+COS2")
        #check = check and (stepsList[1]['repeat'] == 1)
        #self.assertTrue(check)


    def test_getTreadFrameFunction(self):
        '''A very simple case

        outExpected: "LINE101=(0, -3, MASK1)
        '''
        sqNum = "101"
        startPoint = [0, -3]
        stepMask = "MASK1=(o, on, 100, off)"
        slope = 0
        outExpected = "LINE101=(0, -3, MASK1)"
        outGot = functionGenerator.getTreadFunction(sqNum, startPoint, stepMask, slope)
        self.assertEqual(outExpected, outGot)

        '''One with non-zero slope

        outExpected = "LINE101=(0.5, -20, MASK1)"
        '''
        slope = 0.5
        startPoint = [10, -20]
        outExpected = "LINE101=(0.5, -20, MASK1)"
        outGot = functionGenerator.getTreadFunction(sqNum, startPoint, stepMask, slope)
        self.assertEqual(outExpected, outGot)


    def test_getRiserFunction(self):
        '''A very simple case
        outExpected - 'LINE105=(-5, 30, MASK6)'
        '''
        sqNum = '105'
        startPoint = [9, -15]
        connectorMask = 'MASK6=()'
        height = 5
        slope = 0
        length = 1
        outExpected = 'LINE105=(-5.0, 30.0, MASK6)'
        outGot = functionGenerator.getRiserFunction(sqNum, startPoint, length, connectorMask, height, slope)
        self.assertEqual(outExpected, outGot)

        ''' One with non-zero slope
        outExpected = 'LINE106=(-4.5, 30.0, MASK6)'
        '''
        sqNum = '106'
        slope = 0.5
        outExpected = 'LINE106=(-4.5, 30.0, MASK6)'
        outGot = functionGenerator.getRiserFunction(sqNum, startPoint, length, connectorMask, height, slope)
        self.assertEqual(outExpected, outGot)


    def test_convertTreadVariateFunctions(self):
        ''' A very simple case, when no variate functions applied
        outExpected = []
        '''
        seqNum = '30101'
        step_function = 'flat'
        reachPara = {'Length': 50,
                'Slope': 0.5,
                'MASK0': 'ALL',
                'UserFunctions': {'Mask0': 'ALL',
                    'MASK1':['0', 'on', '9', 'off']
                                }
                }
        step_len = 10
        maskFun = 'MASK1=(0, on, 9, off)'
        startPoint = [0, 0]
        outExpected = []
        outGot = functionGenerator.convertTreadVariateFunctions(seqNum, step_function, reachPara, step_len, maskFun, startPoint)
        self.assertEqual(outGot, outExpected)

        ''' A more complicated case, where a cosine function need to be converted.
        outExpected = ['COS301011=(-3.0, 5.0, 0, MASK0)']
        '''
        step_function = 'COS1'
        reachPara['COS1'] = "(-3, 1, 0, 'MASK0')"
        outExpected = ['COS301011=(-3.0, 5.0, 0.0, MASK1)']
        outGot = functionGenerator.convertTreadVariateFunctions(seqNum, step_function, reachPara, step_len, maskFun, startPoint)
        self.assertEqual(outGot, outExpected)


    def test_simulateThalweg(self):
        ''' A very simple case. Staight thalweg, no slope.
        thalwegExpect = [1001] * 10
        '''
        xlen = 10
        slope = 0
        functions = []
        centerExpect = np.array(list(range(10)))
        thalwegExpect = np.array([0]*10)
        centerOut, thalwegOut = functionGenerator.simulateThalweg(centerExpect, xlen, slope, functions)
#        print(centerExpect)
#        print(thalwegOut)
#        print(thalwegExpect)
        self.assertFalse(np.array_equal(centerOut, centerExpect))
        self.assertTrue(np.array_equal(thalwegOut, thalwegExpect))

        ''' A very simple case. Staight thalweg, no slope.
        thalwegExpect = [1001] * 10
        '''
        y_v = np.array([0]*xlen)
        centerExpect = np.array(list(range(10)))
        thalwegExpect = np.array([0]*10)
        centerOut, thalwegOut = functionGenerator.simulateThalweg(y_v, xlen, slope, functions)
#        print(centerOut)
#        print(centerExpect)
#        print(thalwegOut)
#        print(thalwegExpect)
        self.assertTrue(np.array_equal(centerOut, centerExpect))
        self.assertTrue(np.array_equal(thalwegOut, thalwegExpect))

        ''' No additional functions, but slope is non-zero.
        thalwegExpect = [1001, 1000, 999, 998, ..., 992]
        '''
        slope = 1
        ls = list(range(-9, 1))
        ls.reverse()
        thalwegExpect = np.array(ls)
        centerOut, thalwegOut = functionGenerator.simulateThalweg(y_v, xlen, slope, functions)
#        print(thalwegOut)
#        print(thalwegExpect)
        self.assertTrue(np.array_equal(centerOut, centerExpect))
        self.assertTrue(np.array_equal(thalwegOut, thalwegExpect))

        ''' Additional functions, slope is non-zero.
        thalwegExpect = [1001, 999, 997, 995, ..., 982]
        '''
        functions = ['LINE1=(-1, 0, MASK0)']
        ls = [0 - 2*x for x in range(10)]
        thalwegExpect = np.array(ls)
        centerOut, thalwegOut = functionGenerator.simulateThalweg(y_v, xlen, slope, functions)
#        print(thalwegOut)
#        print(thalwegExpect)
        self.assertTrue(np.array_equal(centerOut, centerExpect))
        self.assertTrue(np.array_equal(thalwegOut, thalwegExpect))

        ''' Additional functions, slope is non-zero, arbitrary mask function
        thalwegExpect = [1001, 999, 997, 995, ..., 982]
        '''
        functions = ['LINE1=(-1, 0, MASK1030)']
        ls = [0 - 2*x for x in range(10)]
        thalwegExpect = np.array(ls)
        centerOut, thalwegOut = functionGenerator.simulateThalweg(y_v, xlen, slope, functions)
#        print(thalwegOut)
#        print(thalwegExpect)
        self.assertTrue(np.array_equal(centerOut, centerExpect))
        self.assertTrue(np.array_equal(thalwegOut, thalwegExpect))


    def test_getCenterlineFunctions(self):
        ''' No meandering centerline function in reachPara.
        outExpected=[]
        '''
        reachPara = {}
        outExpected=[]
        outGot = functionGenerator.getCenterlineFunctions(reachPara)
        self.assertListEqual(outGot, outExpected)

        ''' Some meandering centerline functions in reachPara.
        reachPara = {'Meandering Centerline Function':'LINE1+LINE2',
                        'LINE1':'(1, 0, MASK0)',
                        'LINE2':'(2, 0, MASK0)'}
        outExpected = ['LINE1'='(1, 0, MASK0)','LINE2'='(2, 0, MASK0)']
        '''
        reachPara = {'Meandering Centerline Function':'LINE1+LINE2',
                'LINE1':'(1, 0, MASK0)',
                'LINE2':'(2, 0, MASK0)'}
        outExpected = ['LINE1=(1, 0, MASK0)','LINE2=(2, 0, MASK0)']
        outGot = functionGenerator.getCenterlineFunctions(reachPara)
        self.assertListEqual(outGot, outExpected)


    def test_simulateCenterline(self):
        ''' A very simple case. Staight thalweg, no slope.
        centerline expected: [0, 1, ..., 9]
        '''
        xlen = 10
        slope = 0
        functions = []
        sExpected = np.array(list(range(10)))
        yExpected = np.array([0]*xlen)
        slopeExpected = 0
        sOut, yOut, slopeOut = functionGenerator.simulateCenterline(xlen, slope, functions)
#        print(centerOut)
#        print(centerExpect)
        self.assertTrue(np.array_equal(sOut, sExpected))
        self.assertTrue(np.array_equal(yOut, yExpected))
        self.assertEqual(slopeOut, slopeExpected)

        ''' None-zero slope, no additional functions.
        centerline expected: [0, 1, ..., 9]
        '''
        slope = 1
        slopeExpected = 1
        sOut, yOut, slopeOut = functionGenerator.simulateCenterline(xlen, slope, functions)
        self.assertTrue(np.array_equal(sOut, sExpected))
        self.assertTrue(np.array_equal(yOut, yExpected))
        self.assertEqual(slopeOut, slopeExpected)

        ''' None-zero slope, additional functions.
        centerline expected: [0, 1, ..., 9]
        '''
        functions = ['LINE1=(4/3, 0, MASK0)']
        ls = [5/3*x for x in range(10)]
        sExpected = np.array(ls)
        yExpected = np.array([4/3*x for x in range(10)])
        slopeExpected = 0.67
        sOut, yOut, slopeOut = functionGenerator.simulateCenterline(xlen, slope, functions)
        self.assertEqual(np.round(slopeOut, 2), slopeExpected)
        for i in range(len(yExpected)):
            self.assertEqual(np.round(sExpected[i], 4), np.round(sOut[i], 4))
            self.assertEqual(np.round(yExpected[i], 4), np.round(yOut[i], 4))


    def test_calFunction(self):
        ''' A very simple case, no repeat, only one step
        '''
        y_center = np.array([0]*50)
        s_center = np.array(list(range(50)))
        ind = '1'
        step = {'tread length': 9,
                'riser height': 5,
                'riser length': 1,
                'tread function': 'flat',
                'repeat': 1}

        startPoint = [0, -15]
        reachPara = {'Length': 50,
                'riverSlope': 0,
                'MASK0': 1,
                'UserFunctions': {}
                }

        functionsExpected = ['MASK110=(0, on, 9, off)',
                            'LINE110=(0, -15, MASK110)',
                            'MASK112=(9, on, 10, off)',
                            'LINE112=(-5.0, 25.0, MASK112)'
                ]
        startPointExpected = [9, -20.0]

        functions, startPoint = functionGenerator.calFunction(s_center, y_center, ind, step, startPoint, reachPara)
        self.assertEqual(functionsExpected, functions)
        self.assertEqual(startPointExpected, startPoint)

        ''' When a reach repeat more than once.'''
        step['repeat'] = 2
        functionsExpected = ['MASK110=(0, on, 9, off)',
                            'LINE110=(0, -15, MASK110)',
                            'MASK112=(9, on, 10, off)',
                            'LINE112=(-5.0, 25.0, MASK112)',
                            'MASK120=(10, on, 19, off)',
                            'LINE120=(0, -20.0, MASK120)',
                            'MASK122=(19, on, 20, off)',
                            'LINE122=(-5.0, 70.0, MASK122)'
                ]
        startPointExpected = [19.0, -25.0]

        startPoint = [0, -15]
        functions, startPoint = functionGenerator.calFunction(s_center, y_center, ind, step, startPoint, reachPara)
        self.assertEqual(functionsExpected, functions)
        self.assertEqual(startPointExpected, startPoint)

        ''' Error example; step length exceeds.'''
        functionsExpected = []
        startPointExpected = [20.0, -25.0]

        startPoint = [45, -15]
        startPointExpected = [45, -15]
        functions, startPoint = functionGenerator.calFunction(s_center, y_center, ind, step, startPoint, reachPara)
        self.assertEqual(functionsExpected, functions)
        self.assertEqual(startPointExpected, startPoint)

        ''' Error example; connector length exceeds.'''
        functionsExpected = ['MASK110=(40, on, 49, off)',
                            'LINE110=(0, -15, MASK110)'
                            ]

        startPoint = [39, -15]
        startPointExpected = [48, -15.0]
        functions, startPoint = functionGenerator.calFunction(s_center, y_center, ind, step, startPoint, reachPara)
        self.assertEqual(startPointExpected, startPoint)
        self.assertEqual(functionsExpected, functions)


if __name__ == '__main__':
    unittest.main()
