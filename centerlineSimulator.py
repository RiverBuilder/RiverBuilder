'''
This serves for users to test their centerline functions.

User only need to provide the length of the river and functions
of the centerline. The input format file should be in the same folder
of this script called 'centerlineSimulatorInputFormat.txt'.

It will produces a imaged called 'simulatedCenterline' showing the 
resulted centerline from user inputs.
'''

from riverbuilder import river
import sys
import numpy as np
import matplotlib.pyplot as plt

def simulateCenterline(fname):
    '''Simulate a centerline with parameters provided by the user input file.
    
    fname - input file name.

    Output:
    The program will generate an image saved in 'simulatedCenterline' file.
    '''
    paraDict = river.fileParser(fname)
    paraDict, info = river.paraCheck(paraDict, "Length", 1000, "int", 1)
    length = paraDict['Length']
    funDict, log = river.buildFunDict(paraDict)
    info += log+'\n'

    ar = np.array(list(range(length)))
    fun, info = river.buildFun("Centerline Function", paraDict, funDict, river.defaultFunction)
    ar, y = fun(ar)

    fig, ax = plt.subplots(1,1)
    yUpper = max(length/2, np.amax(y))
    yLower = min(-length/2, np.amin(y))
    plt.ylim(yLower, yUpper)
    ax.plot(ar, y)
    plt.savefig('simulatedCenterline.png')

    print(info)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('An input file is needed for this script to run.')
        sys.exit()
    fname = sys.argv[1]
    simulateCenterline(fname)
