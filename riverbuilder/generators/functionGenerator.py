'''
This script serves riverbuilder's user to generate functions
simulating thalweg of a step pool river.

One input file name should be given as parameter.

The input formatting can be found in 'stepPoolFunctionGeneratorFormat.txt'.

The program will output two files:
    steppoolThalwegFunctions.txt -- a txt file containing steppool thalweg functions.
    steppolThalwegOverview.png -- a png image file preview the thalweg.
'''

from ..core import river
from ..core.cChannel import Channel
from . import functionTransformer
import os
import sys
import numpy as np
import matplotlib.pyplot as plt


def reachParasCheck(reachPara):
    '''Validate reach parameters.
    This method will validate the following reach parameters:
        Length - int; default is 100
        Slope - int; default is 0
    
    In:
        reachPara - dict, containing information of reach parameters.
    Out:
        reachPara - dict, containing verified reach parameters.
    '''
    reachPara, dummy = river.paraCheck(reachPara, 'Length', 100, 'int', 1)
    reachPara, dummy = river.paraCheck(reachPara, 'Slope', 0, 'float', 0)
    reachPara, dummy = river.paraCheck(reachPara, 'Hbf', 1, 'float', 1)
    return reachPara


def stepParasCheck(stepPara):
    '''Validate step parameters.
    This method will validate the following step parameters:
            tread length - int, default 10,
            riser height - int, default 2,
            riser length - int, default 1,
            tread function - string, default 'flat',
            repeat - int, default 1
    
    In:
        stepPara - dict, containing information of step parameters.
    Out:
        stepPara - dict, containing verified step parameters.
    '''
    stepPara, dummy = river.paraCheck(stepPara, 'tread length', 10, 'int', 1)
    stepPara, dummy = river.paraCheck(stepPara, 'riser height', 2, 'float', 1)
    stepPara, dummy = river.paraCheck(stepPara, 'riser length', 1, 'int', 1)
    stepPara, dummy = river.paraCheck(stepPara, 'tread function', 'flat', 'str')
    stepPara, dummy = river.paraCheck(stepPara, 'repeat', 1, 'int', 1)

    return stepPara


def parseFile(fname):
    '''
    Read information in fname file, and store it to proper output data.

    Out:
    reachPara - dictionary;
    stepsList - list; elem: dictionary
    '''
    log = ''
    try:
        f = open(fname, "r")
    except IOError:
        log = "Error! "+fname+" is not found in current folder:\n"+os.getcwd()+"\nProgram exits.\n"
        print(log)
        sys.exit()

    reachPara = {}
    stepPara = {}
    stepsList = []
    flagStep = False

    lines = f.readlines()
    for line in lines:
        line = line.strip()
        if line.startswith('#') or line == '':
            continue
        elif line == 'step starts':
            flagStep = True
            stepPara = {}
        elif line == 'step ends':
            flagStep = False
            stepPara = stepParasCheck(stepPara)
            stepsList.append(stepPara)
        elif '=' in line:
            [key, val] = line.split('=')
            if flagStep:
                stepPara[key] = val
            elif key not in reachPara:
                reachPara[key] = val
            else:
                reachPara[key] = reachPara[key] + '+' + val

    reachPara = reachParasCheck(reachPara)
    print("Parsing user defined functions...")
    funcDict, log = river.buildFunDict(reachPara)
    print(log)
    print('')
    reachPara["UserFunctions"] = funcDict
    return reachPara, stepsList


def generateMask(i, start, end):
    '''Return a mask function string.

    In:
    i - string; mask number
    start - mask start point
    end - mask end point

    Out:
    mask - string; 'MASKi=(start, on, end, off)
    '''
    return 'MASK'+i+'=('+str(start)+', on, '+str(end)+', off)'


def getTreadFunction(sqNum, startPoint, treadMask, slope):
    '''Get functions for a flat tread.

    In:
    sqNum - string; sequence number of function
    startPoint - [x, y]; where the tread starts
    treadMask - string; mask function for the tread
    slope - float; slope of the river

    Out:
    functions - String; funtion of a flat tread
    '''
    mask = treadMask.split('=')[0]
    function = "LINE"+sqNum+'=('+str(slope)+', '+str(startPoint[1])+', '+mask+')'
    return function


def convertTreadVariateFunctions(seqNum, tread_function, reachPara, tread_len, maskFun, startPoint):
    '''Convert user defined functions in each tread.

    In:
    seqNum - string; sequence number of function
    tread_function - string; functions user want to apply on this tread. ex: "SIN1+COS2"
    reachPara - dictionary;
                {'Length':int,
                 'Slope': float // this should be river slope
                 'UserFunctions':{}
                 ...
                }
    maskFun - string; mask function for this tread
    startPoint - [int, int]; x, y value for the tread starting point.

    Out:
    functions - list; elem: string, function declaration.
    '''

    if tread_function == 'flat':
        return []

    tread_func = tread_function.split('+')
    reach_len = reachPara['Length']

    functions = []

    for t in range(len(tread_func)):
        funcName = tread_func[t]
        func = river.getfunParas(reachPara[funcName], reachPara['UserFunctions'])
        func = functionTransformer.transform(funcName, func, maskFun, reach_len, tread_len, startPoint)

        # check if the function fails to convert
        if func == '':
            continue

        # give a proper function indexing number.
        func = func.split('=')
        intraInd = intToStrComplete(t+1, len(tread_func))
        func[0] = func[0] + seqNum + intraInd
        func = '='.join(func)

        # append to return
        functions.append(func)

    return functions


def getRiserFunction(sqNum, startPoint, length, riserMask, height, slope):
    '''Get functions for a straight riser.
    
    In:
    sqNum - string; sequence number of function
    startPoint - [x, y]; where the riser starts
    length - int; riser length
    height - float; riser height
    riserMask - string; mask function for the riser
    slope - float; slope of the river

    Out:
    functions - String; funtion of a riser
    '''
    mask = riserMask.split('=')[0]
    riserSlope = -height/length + slope
    offset = startPoint[1] - (riserSlope-slope) * (startPoint[0])
    function = "LINE"+sqNum+'=('+str(riserSlope)+', '+str(offset)+', '+mask+')'
    return function


def indexOf(value, array):
    return np.where(array == value)[0][0]


def calFunction(s_v, y_v, ind, step, startPoint, reachPara):
    '''Basing on the step info, calculate out the functions needed.

    Two things are checked:
        1. Overall step length should be less or equal to reach length.
            -- a violation will ignore the exceeding steps.
        2. A reach always ends with a pool instead of a riser.
            -- a vilation will ignore the last riser.

    In:
    ind - string; index of this step segment
    step - dictionary;
            {'tread length': int,
             'riser height': float,  (riser height is the height of the riser)
             'riser length': int
             'repeat': int
             'tread function': string
            }
    startPoint - [int, float]; x and z value; starting point of this step
    reachPara - dictionary;
                {'Length':int,
                 'Slope': float // this should be river slope
                 'FUN1':[val1, val2, val3, val4]
                 'FUN2':[val1, val2, val3, val4]
                 ...
                }

    Out:
    stepFunctions - list; elem: string
    startPoint - [int, float]; x and z value; ending point of this step
    '''
    length = step['tread length']
    height = step['riser height']
    riser_length = step['riser length']
    repeat = step['repeat']
    user_function = step['tread function']
    stepLength = length + riser_length
    reach_len = reachPara['Length']
    slope = reachPara['riverSlope']

    stepFunctions = []
    indexStart = indexOf(startPoint[0], s_v)
    if indexStart != 0:
        indexStart += 1
    tot_length = indexStart
    doRiser = True
    for i in range(repeat):
        # check if step length exceeds reach
        tot_length += length
        if tot_length >= reach_len:
            print('Accumulated steps length exceeds reach length.')
            return stepFunctions, startPoint

        # check if riser length exceeds reach
        tot_length += riser_length
        if tot_length + 1 > reach_len:
            length = reach_len - indexStart - 1
            doRiser = False

        intraInd = intToStrComplete(i+1, repeat)
        # mask functions
        treadMask = generateMask(ind+intraInd+'0', indexStart, indexStart+length)
        riserMask = generateMask(ind+intraInd+'2', indexStart+length, indexStart+stepLength)

        stepFunctions.append(treadMask)

        # calculate function for a flat step-pool
        frameFunction = getTreadFunction(ind+intraInd+'0', startPoint,  treadMask, slope)
        stepFunctions.append(frameFunction)

        # convert tread variate functions
        variateFunctions = convertTreadVariateFunctions(ind+intraInd+'1', user_function, reachPara, length, treadMask, startPoint)
        stepFunctions += variateFunctions

        # determine end point of the pool
        simulate_x, simulate_y = simulateThalweg(y_v, tot_length-riser_length, slope, [frameFunction]+variateFunctions)
        #startPoint = [int(simulate_x[startPoint[0]+length]), float(simulate_y[startPoint[0]+length])]
        startPoint = [int(simulate_x[-1]), float(simulate_y[-1])]
        indexStart = indexOf(startPoint[0], s_v)

        # calculate riser function
        if doRiser:
            stepFunctions.append(riserMask)
            riserFunction = getRiserFunction(ind+intraInd+'2', startPoint, riser_length, riserMask, height, slope)
            stepFunctions.append(riserFunction)

            # determine end point of the riser
            simulate_x, simulate_y = simulateThalweg(y_v, startPoint[0]+riser_length+1, slope, [riserFunction])
            #startPoint = [int(simulate_x[startPoint[0]+riser_length]), float(simulate_y[startPoint[0]+riser_length])]
            startPoint = [int(simulate_x[-1]), float(simulate_y[-1])]
            indexStart = indexOf(startPoint[0], s_v) +1

    return stepFunctions, startPoint


def simulateCenterline(xlen, slope, functions):
    '''Use reverbuilder to simulate centerline with given reach parameters and functions.

    Parameter:
        xlen - int; reach length
        slope - float; reach slope
        functions - list; elem: function declaration in riverbuilder format

    Output:
        s_v - array; values of centerline points along meandering stream.
        y_v - array; y values of centerline points along x-y coordinate system.
        riverSlope - river slope, not valley slope.
    '''
    log = ''

    # Create function dictionary
    funDict = {}
    funCallDict = {}

    for f in functions:
        [key, value] = f.split('=')
        funDict[key] = value
        if key.startswith('MASK'):
            continue
        if "Callee" not in funCallDict:
            funCallDict["Callee"] = key
        else:
            funCallDict["Callee"] = funCallDict["Callee"] + "+" + key

    funDict, log = river.buildFunDict(funDict)

    # Integrate all called functions together
    macroFun, dummy = river.buildFun("Callee", funCallDict, funDict, river.defaultFunction)

    # Simulate centerline
    channel = Channel(xlen, 0, slope)
    channel.setCenterline(macroFun)

    return channel.getCenterline_sn(), channel.getCenterline_y(), channel.getRiverSlope()


def intToStrComplete(x, xtot):
    '''Return string of x with length equal to length of xtot.
    
    Input:
        x - int
        xtot - int

    Ex:
        1, 3 -> '1'
        1, 33 -> '01'
        1, 123 -> '001'

    Out:
        string
    '''
    digitTot = len(str(xtot))
    digitNum = len(str(x))
    return '0'*(digitTot - digitNum) + str(x)


def simulateThalweg(y_center, xlen, slope, functions):
    '''Use reverbuilder to simulate riverthalweg with given reach parameters and functions.

    Parameter:
        xlen - int; reach length
        slope - float; reach slope
        functions - list; elem: function declaration in riverbuilder format

    Output:
        x, y - numpy arrays representing points in x-z coordinate system.
    '''
    log = ''

    # Create function dictionary
    funDict = {}
    funCallDict = {}

    for f in functions:
        [key, value] = f.split('=')
        funDict[key] = value
        if key.startswith('MASK'):
            continue
        if "Callee" not in funCallDict:
            funCallDict["Callee"] = key
        else:
            funCallDict["Callee"] = funCallDict["Callee"] + "+" + key

    funDict, log = river.buildFunDict(funDict)

    # Integrate all called functions together
    macroFun, dummy = river.buildFun("Callee", funCallDict, funDict, river.defaultFunction)

    # Simulate thalweg
    channel = Channel(xlen, 0, slope)
    channel.setCenterlineManual(y_center)
    channel.setHbfManual(1)
    channel.setThalweg(0, macroFun)

    return channel.getCenterline_sn(), channel.getThalweg()


def intToStrComplete(x, xtot):
    '''Return string of x with length equal to length of xtot.
    
    Input:
        x - int
        xtot - int

    Ex:
        1, 3 -> '1'
        1, 33 -> '01'
        1, 123 -> '001'

    Out:
        string
    '''
    digitTot = len(str(xtot))
    digitNum = len(str(x))
    return '0'*(digitTot - digitNum) + str(x)


def getCenterlineFunctions(reachPara):
    '''
    Return a list of functions to calculate centerline.

    Input:
    reachPara - dict; containing information for centerline functions.
                    {...,
                     'Meandering Centerline Function': 'LINE1+LINE2',
                     'LINE1': '(4/3, 0, MASK0)',
                     ...}

    Output:
    [..., 'LINE1=(4/3, 0, MASK0)', ...]
    '''
    funNames = reachPara.get('Meandering Centerline Function')
    if funNames is None:
        return []
    else:
        funNames = funNames.split('+')
        functions = [ x+'='+reachPara[x] for x in funNames]

    for (key, value) in reachPara.items():
        if key.startswith('MASK'):
            functions.append(key+'='+value)
    return functions


def generateFunctions(fname):
    '''
    Parse an input file containing information of a specific
    step-pool river; calculate can return functions can be used
    to generate thalweg of this step-pool river, and plot info 
    to plot a simulated image of the thalweg.

    Out:
    functions - list; elem: string  
                //each string is a function declaration in riverbuilder
                    inputfile format.
    plot - matplotlib.pyplot plot object
    '''
    reachPara, stepsList = parseFile(fname)
    len_reach = reachPara['Length']
    slope = reachPara['Slope']

    # Simulate centerline
    centerlineFunctions = getCenterlineFunctions(reachPara)
    s_v, y_v, riverSlope = simulateCenterline(len_reach, slope, centerlineFunctions)
    reachPara['riverSlope'] = riverSlope

    # Calculate functions
    functions = []
    startPoint = [0, 0]
    stepCounts = len(stepsList)
    for i in range(stepCounts):
        step = stepsList[i]
        ind = intToStrComplete(i+1, stepCounts)
        stepfunctions, startPoint = calFunction(s_v, y_v, ind, step, startPoint, reachPara)
        functions += stepfunctions

    # Complete the rest of reach
    indexStart = indexOf(startPoint[0], s_v) +1
    if startPoint[0] < len_reach - 1:
        length = len_reach - indexStart
        #length = len_reach - startPoint[0] - 1
        lastMask = generateMask('9999', indexStart, indexStart+length)
        lastFunction = getTreadFunction('9999', startPoint, lastMask, slope)
        functions.append(lastMask)
        functions.append(lastFunction)

    # Throw these functions into riverbuilder and get out a plot.
    x, y = simulateThalweg(y_v, len_reach, slope, functions)
#    x = np.array(list(range(100)))
#    y = np.sin(2*x/99*np.pi)

    fig, ax = plt.subplots(1,1)
    ax.plot(x, y)

    return functions, fig


if __name__ == '__main__':
    '''
    Parse one argument which should be input file.
    Call generateFunctions to get result
    Output results to :
        screen -- generated functions printout.

    '''
    if len(sys.argv) == 1:
        print('An input file is needed for this script to run.')
        sys.exit()
    fname = sys.argv[1]

    functions, fig = generateFunctions(fname)
    functions = [x+'\n' for x in functions]
    outTxt = ["# Copy and Paste the following lines to USER-DEFINED FUNCTIONS section\n\n"]
    outTxt += functions
    outTxt.append("\n")

    outTxt += ["# Copy and Paste the following lines to CHANNEL BREAKLINE INPUT PARAMETERS section\n\n"]
    for function in functions:
        funcName = function.split('=')[0]
        if funcName.startswith("MASK"):
            continue
        outTxt.append("Thalweg Elevation Function="+funcName+'\n')

    try:
        f = open('steppoolThalwegFunctions.txt', 'w')
        f.writelines(outTxt)
    except:
        print('Unable to write to file:', 'steppoolThalwegFunctions.txt\n\
                Instead, print all functions on console.')

    print("The calculated functions are:")
    for fun in functions:
        print(fun[:-1])

    plt.savefig('steppoolThalwegOverview')
