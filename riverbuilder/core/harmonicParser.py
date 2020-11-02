from .functions import *
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import sys


def extractFun(x_v, y_v, funID=0):
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

    out = ','.join([str(amplitude), str(frequency), str(phase), "MASK0"])
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

    print('Input ready')
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
