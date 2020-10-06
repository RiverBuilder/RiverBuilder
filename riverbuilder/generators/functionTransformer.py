'''
This script transform a function from one reach to another.
'''
import re
from math import pi

FUNCNAME = set(["SIN", "COS", "LINE", "PERLIN", "SINSQ", "COSSQ", "CNOIDAL", "STEP"])

def transformTri(function, maskFunc, reach_len, tread_len, startPoint):
    '''Transform Sin function expression.
    
    Input:
        function - list. [val1, val2, val3, val4]
        maskFunc - string; 'MASK1=(1, on, 2, off)'
        reach_len - int;
        tread_len - int;
        startPoint - tuple; (x, y)

    Output:
        func - string; "(amp, frq, shift, Mask)"
    '''
    amplitude = function[0]
    frequency = function[1] * (reach_len/tread_len)
    shift = (function[2] - startPoint[0])*frequency*2*pi/reach_len # this may be wrong
    func = ','.join([str(amplitude), ' '+str(frequency), ' '+str(shift), ' '+maskFunc])
    func = '('+func+')'
    return func


def transformLine(function, maskFunc, reach_len, tread_len, startPoint):
    slope = function[0]
    offset = startPoint[1] - slope*startPoint[0]
    func = ','.join([str(slope), ' '+str(offset), ' '+maskFunc])
    func = '('+func+')'
    return func


def transformPerlin(function, maskFunc, reach_len, tread_len, startPoint):
    func = ''
    return func


def transformCnoidal(function, maskFunc, reach_len, tread_len, startPoint):
    func = ''
    return func


def transformStep(function, maskFunc, reach_len, tread_len, startPoint):
    func = ''
    return func


TRANSFORMER = {'SIN': transformTri,
                'COS': transformTri,
                'LINE': transformLine,
                'PERLIN': transformPerlin,
                'SINSQ': transformTri,
                'COSSQ': transformTri,
                'CNOIDAL': transformCnoidal,
                'STEP': transformStep
                }


def transform(funcName, function, maskFunc, reach_len, tread_len, startPoint):
    '''Transform a user defined function into a proper tread function
    
    Input:
        funcName - string; "SIN3"
        function - list; [1, 1, 0, "MASK1"]
        maskFunc - string; "MASK1=(1, on, 50, off)"
        reach_len - int
        tread_len - int
        startPoint - tuple; (0, 0)

    Output:
        string; "(1, 2.0, 0, MASK2)"
    '''
    # extract fun name
    p = re.compile('[A-Z]+')
    m = p.match(funcName)
    funcName = m.group()
    maskFunc = maskFunc.split('=')[0]

    # validate function name
    if funcName not in FUNCNAME:
        print(funcName,'is not an allowed function function.')
        return ''

    # get the transformed expression of the function
    exp = TRANSFORMER[funcName](function, maskFunc, reach_len, tread_len, startPoint)

    return funcName + '=' + exp
