'''This is the valley class module.

It handles simulating a valley containing a river.
It consists of a river, which is a Channel object, and arbitrary
number of valley levels.

All offsets from valley levels to river are calculated in xy-coordinate system.
'''

import numpy as np
from . import functions
from .cPipe import Pipe
from .cChannel import Channel
from math import pi, sqrt, log
import csv
import matplotlib.pyplot as plt


class Valley(Pipe):

    def __init__(self, channel=Channel()):
        '''Initiator for Valley class

        x_len -- length of valley in x direction
        channel -- Channel object passed to a valley
        x_slope -- slope in x direction
        dx -- resolution in x direction

        Class attributes:
        reshape -- flag show that if the channel has been reshaped
        channelCenter_x -- x values for the points of channel centerline
        channelCenter_y -- y values for the points of channel centerline
        '''
        self.channel = channel
        super().__init__(int(channel.x_len), channel.x_slope, channel.dx, channel.zd)
        self.channelCenter_x = channel.x_v
        self.channelCenter_y = channel.getCenterline_y()
        self.setThalweg(zd=(-1*channel.zd))


    def setValleyBoundary(self, z_offset, z_start, y_offset, direction, yfun=None):
        z_start = np.zeros(len(self.x_v))
        z_max = np.amax(self.channel.levels_z[direction][-1])*self.dx
        z_offset += z_max

        if direction == 'left':
            y_max = np.amax(self.channel.levels_y[direction][-1])*self.dx
            y_offset += y_max
        else:
            y_max = np.amin(self.channel.levels_y[direction][-1])*self.dx
            y_offset = abs(y_max - y_offset)


        super().setLevel(z_offset, z_start, y_offset, direction, yfun)


    def getLevel(self, dictkey, direction, ind=None):
        '''Rewrite the getLevel function. Use the getBound for the first time.'''
        if self.leveldict[dictkey][direction] == []:
            return self.channel.getLevel(dictkey, direction, ind)
        else:
            return super().getLevel(dictkey, direction, ind)


    def tolist_levelxy(self):
        '''Return x, y values for all levels in a secondary list.'''
        x = [] + self.channelCenter_x.tolist()
        y = [] + self.channelCenter_y.tolist()

        self.helpAppend(x, self.channel.levels_x)
        self.helpAppend(y, self.channel.levels_y)
        self.helpAppend(x, self.levels_x)
        self.helpAppend(y, self.levels_y)

        return [x, y]


    def tolist_levelxz(self):
        '''Return x, y values for all levels in a secondary list.'''
        x = [] + self.channelCenter_x.tolist()
        z = [] + self.channel.getThalweg().tolist()

        self.helpAppend(x, self.channel.levels_x)
        self.helpAppend(z, self.channel.levels_z)
        self.helpAppend(x, self.levels_x)
        self.helpAppend(z, self.levels_z)

        return [x, z]


    def tolist_all(self):
        '''Return x, y z values for everythin in a secondary list.'''

        xshape_x, xshape_y, xshape_z = self.channel.getXShape()
        center_z = self.channel.getCenterlineElevation()
        x = [] + self.channelCenter_x.tolist()
        y = [] + self.channelCenter_y.tolist()
        z = [] + center_z.tolist()
        label = ['CL'] * len(x)

        x += xshape_x.tolist()
        y += xshape_y.tolist()
        z += xshape_z.tolist()
        label += ['XS']*len(xshape_x)

        self.helpAppend(x, self.channel.levels_x)
        self.helpAppend(y, self.channel.levels_y)
        self.helpAppend(z, self.channel.levels_z)
        label += ['RB'] * (len(x) - len(label))

        self.helpAppend(x, self.levels_x)
        self.helpAppend(y, self.levels_y)
        self.helpAppend(z, self.levels_z)
        label += ['VL'] * (len(x) - len(label))


        return [x, y, z, label]


    def getXShapePlot(self):
        '''return matplotlib plot object that contains X-Shape plots of the valley.'''
        slope = self.getSlope()
        minInd = np.argmin(np.absolute(slope))
        minInd_x = self.x_v[minInd]
        
        y = []
        z = []
        for i in range(0, len(self.levels_y['left'])):
            y.append(self.levels_y['left'][-1*i-1][minInd]*self.dx)
            z.append(self.levels_z['left'][-1*i-1][minInd]*self.dx)

        for i in range(0, len(self.channel.levels_y['left'])):
            ind = np.argmin(np.absolute(self.channel.levels_x['left'][-1*i-1] - minInd_x))
            y.append(self.channel.levels_y['left'][-1*i-1][ind]*self.dx)
            z.append(self.channel.levels_z['left'][-1*i-1][ind]*self.dx)

        indLeft = np.argmin(np.absolute(self.channel.levels_x['left'][0] - minInd_x))
        indRight = np.argmin(np.absolute(self.channel.levels_x['right'][0] - minInd_x))
        wbf = self.channel.levels_y['left'][0][indLeft] - self.channel.levels_y['right'][0][indRight]
        if self.channel.tz == -1:
            channel_y, channel_z = self.channel.pointXShape(ind, 0, wbf, self.channel.xshapePoints)
        else:
            channel_y, channel_z = self.channel.suXShape(ind, wbf, self.channel.tz, self.channel.xshapePoints)
        channel_y = channel_y + (self.channel.levels_y['left'][0][indLeft] + self.channel.levels_y['right'][0][indRight])/2

        channel_y = channel_y*self.dx
        channel_z = channel_z*self.dx

        y += channel_y.tolist()
        z += channel_z.tolist()

        for i in range(0, len(self.channel.levels_y['right'])):
            ind = np.argmin(np.absolute(self.channel.levels_x['right'][i] - minInd_x))
            y.append(self.channel.levels_y['right'][i][ind]*self.dx)
            z.append(self.channel.levels_z['right'][i][ind]*self.dx)

        for i in range(0, len(self.levels_y['right'])):
            y.append(self.levels_y['right'][i][minInd]*self.dx)
            z.append(self.levels_z['right'][i][minInd]*self.dx)

        fig, ax = plt.subplots(1, 1)
        fig.suptitle('Valley X-Shape')
        ax.plot(y, z, '-', marker='o')
        plt.xlabel('Y')
        plt.xlabel('Z')
        return fig


    def tocsv(self, outfile):
        '''Outwrite xy values and xz values of all levels to output file.'''
        header = ["X", "Y", "Z", "Label"]
        out = [header]
        xyz = self.tolist_all()

        xyz_out = [[round(xyz[0][i]*self.dx,3), round(xyz[1][i]*self.dx, 3), round(xyz[2][i]*self.dx, 3), xyz[3][i]] for i in range(len(xyz[0]))]
#        xyz_out = [[round(xyz[0][i],3), round(xyz[1][i], 3), round(xyz[2][i], 3), xyz[3][i]] for i in range(len(xyz[0]))]
        out += xyz_out

        with open(outfile+".csv", 'w') as cf:
            cw = csv.writer(cf, lineterminator='\n')
            cw.writerows(out)


    def getValleySlope(self):
        '''Return Valley slope.'''
        return self.getPipeSlope()


    def __str__(self):
        '''Turn most important data of valley.'''
        sl = self.getSL()
        slope = self.getValleySlope()
        s = 'Sinuosity:'+str(round(sl, 3))+'\n'
        s += 'Valley Slope:'+str(slope)+'\n'

        for i in range(0, len(self.levels_n['left'])):
            s += 'Average Width of '+'L'+str(i)+' Valley Level is: '+str(round(np.average(self.levels_n['left'][i])*self.dx,3)) + '\n'
            
        for i in range(0, len(self.levels_n['right'])):
            s += 'Average Width of '+'R'+str(i)+' Valley Level is: '+str(abs(round(np.average(self.levels_n['right'][i]*self.dx), 3))) + '\n'

        return s
###########################################################

    def helpAppend(self, li, dic):
        for array in dic["left"]:
            li += array.tolist()
        for array in dic["right"]:
            li += array.tolist()
