'''This is the Channel module

It can simulate a river channel basing on the inputs it is provided.
It consists of a centerline, an inner channel, and arbitrary number of
outer banks.
All functions apply to it should be continuous.
The offsets from banks to centerline are in sn coordinate system, and 
transform into xy coordinate system later on.
'''

import numpy as np
from . import functions
import math
from math import pi, sqrt, log, ceil, floor
import csv
from .cPipe import Pipe
import matplotlib.pyplot as plt

class Channel(Pipe):

    def __init__(self, x_len=100, wbf_min=0, valley_slope=0.01, dx=1):
        '''Channel class initiator

        x_len -- int; valley length in x direction
        wbf_min -- float; minimum bankfull width
        valley_slope -- float; slope of valley
        dx -- int; resolution in x direction

        class private variables:
        hbf -- float; average bankfull height
        thalweg -- array; z values of thalweg
        curvature -- array; curvature of centerline
        xshapePoints -- int; number of points in each Xshape
        xshape_x -- array; x values of Xshape
        xshape_y -- array; y values of Xshape
        xshape_z -- array; z values of Xshape
        z_center -- array; z values of centerline
        dynamicCurv -- array; values of curvature of center line
        tz -- int; trapezoid xshape bottom points. -1 means asymetric
        '''
        super().__init__(int(x_len), valley_slope, int(dx))
        self.wbf_min = wbf_min
        self.turns_center = []
        self.hbf = None
        self.thalweg = None
        self.curvature = None
        self.xshapePoints = 21
        self.xshape_x = None
        self.xshape_y = None
        self.xshape_z = None
        self.z_center = None
        self.dynamicCurv = None
        self.tz = -1


    def shapeCenterline(self, fun):
        '''Shape the centerline. Basically recalculate the centerline.'''
        x_v = self.x_v
        n_v = self.getCenterline_y()

        x_v_valley = list(set(x_v.tolist()))
        x_v_valley.sort()
        x_v_valley = np.array(x_v_valley)
        x_v_valley, y_v = fun(x_v_valley)
        x_max = np.amax(x_v_valley)

        out_x, out_y = [], []
        for i in range(len(x_v)):
            x = x_v[i]
            ind = np.where(x_v_valley == x)[0][0]
            if x == x_max:
                continue
            x1, x2 = x_v_valley[ind], x_v_valley[ind+1]
            y1, y2 = y_v[ind], y_v[ind+1]
            x_new, y_new = functions.sn_to_xy(x1, y1, x2, y2, n_v[i])
            out_x.append(x_new)
            out_y.append(y_new)
        self.x_v = np.array(out_x)
        self.y_center = np.array(out_y)


    def getRiverSlope(self):
        '''Return river slope'''
        return self.getPipeSlope()


    def setXShapePoints(self, n):
        '''Set how many points in one x-section shape.'''
        self.xshapePoints = n


    def setHbfManual(self, hbf):
        '''Mannually set self.hbf'''
        self.hbf = hbf


    def setHbf(self, d50=0.01, css=0.047, g_s=1922, g_w=1000):
        '''Automatically calculate Hbf'''
        self.hbf = functions.shields(d50, css, self.getRiverSlope(), g_s, g_w)


    def getHbf(self):
        '''Return self.hbf'''
        if self.hbf is None:
            self.setHbf()
        return self.hbf


    def setTZ(self, n):
        self.tz = n


    def setThalweg(self, zd=1000, fun=np.vectorize(lambda x : 0)):
        '''Calculate z values for thalweg

        zd -- Datum
        fun -- function used to calculate thalweg

        Value will be modified:
        self.thalweg
        '''
        hbf = self.getHbf()
        s_v = self.getCenterline_sn()
        x, y = fun(s_v)
        self.thalweg = y - self.getRiverSlope()*s_v + zd
#        self.thalweg = hbf*y + hbf - self.getRiverSlope()*s_v + zd


    def getThalweg(self):
        '''Return self.thalweg'''
        if self.thalweg is None:
            self.setThalweg()
        return self.thalweg


    def setCurvature(self, fun):
        '''Mannually set centerline curvature

        fun -- function to calculate curvature
        '''
        x = self.x_v
        dummy, self.dynamicCurv = fun(x)


    def getDynamicCurv(self):
        '''Return self.dynamicCur'''
        if self.dynamicCurv is None:
            self.setDynamicCurv()
        return self.dynamicCurv


    def createInnerChannel(self, leftFun=None, rightFun=None, zd=1000, thalwegFun=np.vectorize(lambda x : 0)):
        '''Create most inner channel of river

        leftFun -- function to calculate left inner bank
        rightFun -- function to calculate right innert bank
        zd -- Datum
        thalwegFun -- function to calculate thalweg

        Value will be modified:
        self.levels_x
        self.levels_y
        self.levels_z
        self.levels_n
        '''
        self.setThalweg(zd, thalwegFun)

        thalweg = self.getThalweg()
        thalweg_max = np.amax(thalweg)
        flat_thalweg = thalweg + self.getRiverSlope()*self.s_center
        diff = np.amax(flat_thalweg) - flat_thalweg[0]
        z_start = thalweg_max + diff
        hbf = self.getHbf()
        wbf = self.wbf_min/2

        self.setLevel(hbf, z_start, wbf, 'left', leftFun, True)
        self.setLevel(hbf, z_start, wbf, 'right', rightFun, True)


    def getHeight(self, direction):
        '''Return z value of most outer bank of river

        direction -- "left" or "right"
        '''
        if self.levels_z[direction] == []:
            self.createInnerChannel()
        return np.amax(self.levels_z[direction][-1])


    def getAveWbf(self):
        '''Return average bankfull width.'''
        if self.levels_y['left'] == []:
            self.createInnerChannel()
        bf = self.levels_n['left'][0] + np.absolute(self.levels_n['right'][0])
        return np.average(bf)


    def getAveHbf(self):
        '''Return average bankfull height.'''
        if self.levels_y['left'] == []:
            self.createInnerChannel()
        thalweg = self.getThalweg()
        thalweg_max = np.amax(thalweg)
        return np.average(thalweg_max-thalweg) + self.getHbf()


    def getCoWbf(self):
        '''Return coefficient of variation of bankfull width.'''
        ave = self.getAveWbf()
        std = np.std(self.levels_n['left'][0] + np.absolute(self.levels_n['right'][0]))
        return std/ave


    def getCoHbf(self):
        '''Return coefficient of variation of bankfull width.'''
        ave = self.getAveHbf()
        thalweg = self.getThalweg()
        thalweg_max = np.amax(thalweg)
        std = np.std(thalweg_max-thalweg+self.getHbf())
        return std/ave


    def getXShape(self):
        '''Return x, y, z values for Xshape of the whole channel'''
        if self.xshape_x is None:
            self.setXShape()
        return self.xshape_x, self.xshape_y, self.xshape_z


    def getCenterlineElevation(self):
        '''Return z values for centerline.'''
        if self.xshape_x is None:
           self.setXShape()
        return self.z_center


    def getXShapePlot(self):
        '''return matplotlib plot object that contains X-Shape plots of the river Channel.'''
        cur_v = self.getDynamicCurv()
        maxCur = np.amax(cur_v)
        minCur = np.amin(cur_v)

        # If no curvature at all, plot at middle point
        if maxCur == minCur or self.tz != -1:
            fig, ax = plt.subplots(1, 1)
            fig.suptitle('X-Shape for Channel')
            midInd = floor(len(self.x_v)/2)
            wbf = abs(self.levels_n["left"][0][midInd]) + abs(self.levels_n["right"][0][midInd])
            if self.tz == -1:
                y, z = self.pointXShape(midInd, maxCur, wbf, self.xshapePoints)
            else:
                y, z = self.suXShape(midInd, wbf, self.tz, self.xshapePoints)
            z = z + midInd*self.x_slope
            y, z = self.addBankPoints(y, z, midInd)
            ax.plot(y, z, 'k-', marker='o', label='x = '+str(midInd))
            plt.xlabel('Y (related to center of channel)')
            plt.ylabel('Z')
            plt.legend()
            return fig
        else:
            abs_cur_v = np.absolute(cur_v)
            fig, ax = plt.subplots(2, 1, sharex=True)
            fig.suptitle('Max Curvature X-Shape vs. Zero Curvature X-Shape')

            plt.subplot(212)
            indMax = np.argmax(abs_cur_v)
            maxCur = cur_v[indMax]
            wbf = abs(self.levels_n["left"][0][indMax]) + abs(self.levels_n["right"][0][indMax])
            y, z = self.pointXShape(indMax, maxCur, wbf, self.xshapePoints)
            si = self.getCenterline_sn()[indMax]
            z = z + si*self.getPipeSlope()
            y, z = self.addBankPoints(y, z, indMax)
            plt.plot(y, z, 'k-', marker='o', label='Max Curvature:\nx = '+str(indMax))
            plt.xlabel('Y (related to center of channel)')
            plt.ylabel('Z')
            plt.legend()

            plt.subplot(211)
            indMin = np.argmin(abs_cur_v)
            minCur = cur_v[indMin]
            wbf = abs(self.levels_n["left"][0][indMin]) + abs(self.levels_n["right"][0][indMin])
            y, z = self.pointXShape(indMin, maxCur, wbf, self.xshapePoints)
            si = self.getCenterline_sn()[indMin]
            z = z + si*self.getPipeSlope()
            y, z = self.addBankPoints(y, z, indMin)
            plt.plot(y, z, 'k-', marker='o', label='Min Curvature:\nx = '+str(indMin))
            plt.ylabel('Z')
            plt.legend()
            return fig


    def setXShape(self, n=-1):
        '''Calculate x, y, z values for Xshape of the whole channel.
           Also calculate the z values of centerline.
           xshapePointsDict: {(x, y): [z, (x_center, y_center)]}
        '''
        out_x = np.array([])
        out_y = np.array([])
        out_z = np.array([])
        xshapePointsDict = {}
        center_z = []
        y_center = self.getCenterline_y()
        s_center = self.getCenterline_sn()
        pipe_slope = self.getPipeSlope()

        cur_v = self.getDynamicCurv()
        maxCur = np.amax(np.absolute(cur_v))
        for ind in range(len(y_center)-1):
            wbf = abs(self.levels_n['left'][0][ind]) + abs(self.levels_n['right'][0][ind])
            centerOffset = (self.levels_n['left'][0][ind] + self.levels_n['right'][0][ind])/2
            x1 = self.x_v[ind]
            y1 = y_center[ind]

            x2 = self.x_v[ind+1]
            y2 = y_center[ind+1]

            s = s_center[ind]
            # This if statement will determine whether it is AU or SU
            if n == -1:
                y_temp, z = self.pointXShape(ind, maxCur, wbf, self.xshapePoints)
            else:
                y_temp, z = self.suXShape(ind, wbf, n, self.xshapePoints)

            y_temp = y_temp + centerOffset
            real_x, real_y = functions.sn_to_xy(x1, y1, x2, y2, y_temp)
            z = z - pipe_slope*s                                            # use s instead of x

            xshapePointsDict = self.updateXShapePointsDict(xshapePointsDict, real_x, real_y, z, x1, y1)

            #find z for center line
            center_z.append(self.calCenter_z(real_x, real_y, z, x1, y1))

        center_z.append(center_z[-1])

        out_x, out_y, out_z = [], [], []
        for (x, y) in xshapePointsDict.keys():
            out_x.append(x)
            out_y.append(y)
            out_z.append(xshapePointsDict[(x, y)][0])

        out_x = np.array(out_x)
        out_y = np.array(out_y)
        out_z = np.array(out_z)

        self.xshape_x = out_x
        self.xshape_y = out_y
        self.xshape_z = out_z
        self.z_center = np.array(center_z)


    def tolist(self):
        '''Return x, y, z values for all levels in a secondary list'''
        x = []
        y = []
        z = []

        x += self.x_v.tolist()
        y += self.getCenterline_y().tolist()
        z += self.getThalweg().tolist()

        self.helpAppend(x, self.levels_x)
        self.helpAppend(y, self.levels_y)
        self.helpAppend(z, self.levels_z)

        return [x, y, z]


    def tocsv(self, outfile):
        '''Outwrite xy values and xz values of all banks to output file.'''
        header = ["x", "y", 'z']
        out = [header]
        xyz = self.tolist()
        xyz_out = [[round(xyz[0][i], 3), round(xyz[1][i], 3), round(xyz[2][i], 3)] for i in range(len(xyz[0]))]
        out += xyz_out
        with open(outfile+".csv", 'w') as cf:
            cw = csv.writer(cf)
            cw.writerows(out)


    def __str__(self):
        '''Turn most important data of river.'''
        sl = self.getSL()
        aveWbf = self.getAveWbf()
        aveHbf = self.getAveHbf()
        slope = self.getRiverSlope()
        coWbf = self.getCoWbf()
        coHbf = self.getCoHbf()
        s = 'Sinuosity:'+str(round(sl, 3))+'\n'
        s += 'Channel Slope:'+str(slope)+'\n'
        s += 'Average Width of Inner Channel:'+str(round(aveWbf, 3))+'\n'
        s += 'Average Height of Inner Channel:'+str(round(aveHbf, 3))+'\n'
        s += 'Coefficient of Variation (W_ic):'+str(round(coWbf, 3))+'\n'
        s += 'Coefficient of Variation (H_ic):'+str(round(coHbf, 3))+'\n'

        if len(self.levels_n['left']) == 1:
            return s
        for i in range(1, len(self.levels_n['left'])):
            s += 'Average Width Offset of '+'L'+str(i)+' Outer Bank is: '+str(round(np.average(self.levels_n['left'][i]), 3)) + '\n'
            
        if len(self.levels_n['right']) == 1:
            return s
        for i in range(1, len(self.levels_n['right'])):
            s += 'Average Width Offset of '+'R'+str(i)+' Outer Bank is: '+str(abs(round(np.average(self.levels_n['right'][i]), 3))) + '\n'

        return s


##############################################################

    def helpAppend(self, li, dic):
        for array in dic["left"]:
            li += array.tolist()
        for array in dic["right"]:
            li += array.tolist()


    def pointXShape(self, ind, maxCur, wbf, n):
        '''Return y values and z values of XSection of given x

        n -- number of points to calculate XSection
        '''
        cur = self.getDynamicCurv()[ind]
        pipe_slope = self.getPipeSlope()
        si = self.getCenterline_sn()[ind]

        if maxCur == 0:
            B = 1/2
        else:
            B = 1/2 * (1 - abs(cur/maxCur))

        if B < 0.1:
            B = 0.1

        if B == 1:
            L = 1
        else:
            L = -1*log(2)/log(B)

        lbx = self.levels_x['left'][0]
        lb = self.levels_z['left'][0]
        rbx = self.levels_x['right'][0]
        rb = self.levels_z['right'][0]
        lb_ind = np.where(lbx == ind)[0]
        if len(lb_ind) == 0:
            bankH = rb[np.where(rbx == ind)[0][0]]
        else:
            bankH = lb[lb_ind[0]]
        hbf = bankH - self.thalweg[ind]

        n_y = np.array([-wbf*x/(n+1) + (1/2)*wbf  for x in range(1, n+1)])
        Y = (wbf/2-n_y) / wbf

        if cur <= 0:
            n_z = 4 * hbf * (Y**L) * (1-Y**L)
        else:
            n_z = 4 * hbf * ((1-Y)**L) * (1-(1-Y)**L)

        n_z = self.thalweg[ind] + hbf - n_z
        return n_y, n_z


    def suXShape(self, ind, wbf, tzn, n):
        '''Return y values and z values of Symmetric XSection of given x.
        
        ind - index of x value on centerline
        tzn - number of points on base
        n - number of points on XS
        '''
        lbx = self.levels_x['left'][0]
        lb = self.levels_z['left'][0]
        rbx = self.levels_x['right'][0]
        rb = self.levels_z['right'][0]
        lb_ind = np.where(lbx == ind)[0]
        if len(lb_ind) == 0:
            bankH = rb[np.where(rbx == ind)[0][0]]
        else:
            bankH = lb[lb_ind[0]]
        hbf = bankH - self.thalweg[ind]
        n_y = np.array([-wbf*x/(n+1) + (1/2)*wbf  for x in range(1, n+1)])

        n_z = []
        sidePoints = floor(((n-tzn)/2) + 1)

        if (n-tzn) % 2 != 0:
            tzn += 1

        for i in range(1, sidePoints):
            n_z.append((i/sidePoints)*hbf)

        zEnd = n_z[::-1]
        n_z += [hbf] * tzn
        n_z += zEnd

        n_z = self.thalweg[ind] + hbf - n_z
        return n_y, np.array(n_z)


    def addBankPoints(self, y, z, ind):
        '''Add bank points to xshape points.
        
        y - y values for xshape points
        z - z values for xshape points
        ind - where the xshape points are calculated

        Return:
        y, z with bank points added
        '''
        leftEdge = y[0]-(y[1]-y[0])
        y = np.append(leftEdge, y)
        rightEdge = y[-1]+(y[1]-y[0])
        y = np.append(y, rightEdge)

        z = np.append(self.levels_z['left'][0][0], z)
        z = np.append(z, self.levels_z['left'][0][0])

        for i in range(1, len(self.levels_n['left'])):
            y = np.append(self.levels_n['left'][i][ind] - self.levels_n['left'][0][ind] + leftEdge, y)
            z = np.append(self.levels_z['left'][i][0], z)


        for i in range(1, len(self.levels_n['right'])):
            y = np.append(y, self.levels_n['right'][i][ind] - self.levels_n['right'][0][ind] + rightEdge)
            z = np.append(z, self.levels_z['right'][i][0])
        return y, z


    def setDynamicCurv(self):
        ''' Calculate the dynamic curve of the centerline.
        '''
        x_v = self.getCenterline_x()
        y_v = self.getCenterline_y()
        slopeVectorList = []

        for i in range(len(x_v)-1):
            v = (x_v[i+1]-x_v[i], y_v[i+1]-y_v[i])
            slopeVectorList.append(v)

        cur = []
        piCheck = pi/2
        for i in range(len(slopeVectorList)-1):
            v1 = slopeVectorList[i]
            v2 = slopeVectorList[i+1]
            angle = functions.angle_between(v1, v2)
            if np.cross(v1, v2) >= 0:
                cur.append(functions.angle_between(v1, v2))
            else:
                cur.append(functions.angle_between(v1, v2) * -1)

        cur.append(cur[-1])
        cur.append(cur[-1])
        self.dynamicCurv = np.array(cur)


    def calCenter_z(self, real_x, real_y, z, x1, y1):
        '''Calculate the z value for the centerline point.'''
        xpoints = [(real_x[i], real_y[i]) for i in range(len(real_x))]
        minInd = 0
        minDist = functions.pointDist(xpoints[0], (x1, y1))
        diff = [minDist]
        for i in range(1, len(xpoints)):
            dist = functions.pointDist(xpoints[i], (x1, y1))
            diff.append(dist)
            if dist < minDist:
                minInd = i
                minDist = dist

        if minDist == 0 or (minInd-1 < 0 and minInd+1 >= len(diff)) :
            return z[minInd]
        elif minInd+1 < len(diff) and (minInd-1 < 0 or diff[minInd-1] >= diff[minInd+1]):
            minInd2 = minInd+1
        else:
            minInd2 = minInd-1
        
        z1 = z[minInd]
        z2 = z[minInd2]
        minX1 = real_x[minInd]
        minX2 = real_x[minInd2]
        if minX1 == minX2:
            return (z1+z2)/2
        elif minX1 < minX2:
            alpha = (x1-minX1)/(minX2-minX1)
            return alpha*z1 + (1-alpha)*z2
        else:
            alpha = (x1-minX2)/(minX1-minX2)
            return alpha*z1 + (1-alpha)*z2

    
    def updateXShapePointsDict(self, pointsDict, x_v, y_v, z_v, x_center, y_center):
        '''Update the XShape points in XShape points Dictionary.'''
        for i in range(len(x_v)):
            p = (int(x_v[i]), int(y_v[i]))

            if p not in pointsDict:
                pointsDict[p] = [z_v[i], (x_center, y_center)]
                continue

            center_pre = pointsDict[p][1]
            center_current = (x_center, y_center)

            dist_pre = functions.pointDist(p, center_pre)
            dist_current = functions.pointDist(p, center_current)

            if dist_pre > dist_current:
                pointsDict[p] = [z_v[i], (x_center, y_center)]

        return pointsDict