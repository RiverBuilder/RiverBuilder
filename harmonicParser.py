'''
Takes 2 arguments;
    1. Input csv file containing harmonic function data.
        Header should look like: "x, raw_series, harmonic_y1, harmonic_y2,..., all_x_harmonics"
    2. Output file name.

The program will output two files. One is the txt file contains all function strings. The other 
    is a .png file shows the simulated curve.

An example command is like:
    python harmonicParser.py W_base_ft_harmonics_by_power.csv W_base.txt
'''


import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import sys


def slope(x1, y1, x2, y2):
    '''Return slope of two points.'''
    if x1 != x2:
        return (y2-y1) / (x2-x1)
    elif y2 > y1:
        return math.inf
    elif y2 < y1:
        return -math.inf
    else:
        return 0


def slope_v(x_v, y_v):
    '''Return an array of slope.
        si = (y(i+1) - y(i))/(x(i+1) - x(i))
        s(-1) will be the same as s(-2) to make it the same length.
    
    Inputs:
    x_v - array. x values.
    y_v - array. y values.

    Output:
    s_v - array. slopes.
    '''
    x1_v = x_v[:-1]
    x2_v = x_v[1:]

    y1_v = y_v[:-1]
    y2_v = y_v[1:]

    fun = np.vectorize(slope)
    s_v = fun(x1_v, y1_v, x2_v, y2_v) 
    s_v = np.append(s_v, s_v[-1:])
    return s_v


def modx(array):
    '''Return an array is in range [0, 2*pi]'''
    xmin = np.amin(array)
    xmax = np.amax(array)
    xlen = xmax-xmin
    return 2*math.pi*array/xlen

def sin_v(array, amplitude=1, frequency=1, shift=0):
    '''Return the sin values of the array.

    Input:
    array - the values in this array will be transform into 2*pi*xi/max_x
            before calculation.

    Output:
    array - the original array
    out_y - an array of sin values of the original array.
    '''
    array_mod = modx(array)*frequency
    out_y = amplitude*np.sin(array_mod+shift)
    return array, out_y


def extractFun(x_v, y_v, funID=0):
    '''
    return 
    out - str; "SIN[funID]=(amplitude, frequency, phase, "MASK0")"
                which mimic the curve represented by x_v, y_v
    '''
    slope = slope_v(x_v, y_v)
    slope_pre = slope[0:-1]
    slope_post = slope[1:]
    slope_change = np.multiply(slope_pre, slope_post)

    frequency = len(np.where(slope_change <= 0)[0])/2

    y_min = np.amin(y_v)
    y_max = np.amax(y_v)
    if y_min == y_max:
        return '', np.zeros(len(x_v))

    amplitude = (y_max - y_min)/2

    diff_pre = math.inf
    phase = 0
    inc = 0.01
    count = 0
    cp = 0
    while True:
        dummy, simulate = sin_v(x_v, amplitude, frequency, phase)
        diff = np.subtract(y_v, simulate)
        diff = np.sum(np.absolute(diff))

        if diff_pre < diff:
            inc = -inc
            if count > 3 and count - cp <= 3:
                phase -= inc
                break
            cp = count

        phase += inc
        diff_pre = diff
        count += 1

    out = ','.join([str(amplitude), str(frequency), str(round(phase, 4)), "MASK0"])
    out = 'SIN'+str(funID)+'=('+ out+')'

    out_x, out_y = sin_v(x_v, amplitude, frequency, phase)
    return out, out_y


if __name__ == '__main__':
    # Read data
    x_v, y_v_multi = [], []
    with open(sys.argv[1], 'r') as cf:
        rd = csv.reader(cf)
        hd = next(rd)
        y_v_multi = [[] for i in range(len(hd)-3)]
        for row in rd:
            if len(row) == 0 or row[0].strip() == '':
                continue
            x = float(row[0])
            x_v.append(x)

            for i in range(2, len(row)-1):
                y = float(row[i])
                y_v_multi[i-2].append(y)

    # extract sin function
    out = []
    x_v = np.array(x_v)
    y_acc = np.zeros(len(x_v))
    count = 0
    for y_v in y_v_multi:
        y_v = np.array(y_v)
        fun, y_v = extractFun(x_v, y_v, count)

#        fig, ax = plt.subplots(1, 1)
#        ax.plot(x_v, y_v, 'k-')
#        plt.savefig('./harmonic'+str(count))
#
        count +=1
        y_acc = y_acc + y_v
        print(fun)
        out.append(fun+'\n')

    # output txt file
    try:
        f = open(sys.argv[2], 'w')
        f.writelines(out)
        f.close()
    except:
        print('cannot write to ', sys.argv[2])

    fig, ax = plt.subplots(1, 1)
    ax.plot(x_v, y_acc, 'k-')
    plt.savefig('./harmonic')
