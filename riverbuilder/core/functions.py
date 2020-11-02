'''This module contains functions that will be used in riverbuilder project.'''

import numpy as np
import math
from random import random
import random as rd
import copy

def modx(array):
    '''Return an array is in range [0, 2*pi]'''
    xmin = np.amin(array)
    xmax = np.amax(array)
    xlen = xmax-xmin
    return 2*math.pi*array/xlen

# linear function
def line_v(array, slope, intercept, mask):
    '''Return an array of y values based on a linear calculation.'''
    mask = maskCheck(array, mask)
    out = slope*array+intercept
    dummy, out = maskCalculate(out, array, mask)
    return array, out


# mask check
def maskCheck(array, mask):
    '''Check if a mask is valid.
    
    A valid mask should be in form: [DIST1, 0/1, DIST2, 0/1, ...];
    The max dist should be smaller or equal to max value in array.
    The DIST should be increasing in each declaration.
    1 means on in the corresponding segment; 0 means off.
    A special mask is [ALL] indicates always on.

    Output:
    A valid mask with form:
    [(len1, True/False), (len2, True/False), ...]
    '''
    arr_len = len(array)
    default_out = [[0, arr_len, True]]

    if len(mask) == 1:
        return default_out

    xmax = np.amax(array)
    out = []

    #check if mask start with 0
    if mask[0][0].strip() != '0':
        if mask[1].strip().lower() == 'off':
            out.append([0, None, True])
        else:
            out.append([0, None, False])

    for i in range(len(mask)):
        #check distance parameter
        if i % 2 == 0:
            dist = mask[i].strip()
            if not dist.isdigit():
                return default_out

            dist = int(dist)
            #dist = int(int(dist)/xmax*arr_len)
            if len(out) != 0:
                if dist <= out[-1][0]:
                    i += 1  #skip the following switch
                    continue

                if dist >= arr_len:
                    break

                out[-1][1] = dist
                out.append([dist, None, None])
            else:
                out.append([0, None, None])

        # check switch parameter
        else:
            on = mask[i].strip().lower()

            if on == 'off':
                out[-1][-1] = False
            else:
                out[-1][-1] = True


    out[-1][1] = arr_len
    return out


def maskCalculate(y_v, x_v, mask):
    '''Turn on and off on function calculation.

    Inputs:
    y_v - array. y values.
    x_v - array, x values.
    mask - list. [[len1, len2, True/False], [len2, len3, True/False],...]

    Return:
    an undated array with y values.
    '''
    out_y = []
    out_x = []
    default_y = 0

    if len(mask) == 1:
        return x_v, y_v

    for rang in mask:
        start = rang[0]
        end = rang[1]
        on = rang[2]

        #print(start, end)
        for i in range(start, end):
            if on:
                out_x.append(x_v[i])
                out_y.append(y_v[i]+default_y)
            else:
                # reset x:
                a = (i-start) / (end-start)
                if end == len(x_v):
                    a = (i-start) / (end-start-1)
                    out_x.append(x_v[start]*(1-a) + x_v[end-1]*a)
                else:
                    a = (i-start) / (end-start)
                    out_x.append(x_v[start]*(1-a) + x_v[end]*a)
                out_y.append(default_y)

    return np.array(out_x), np.array(out_y)


def sin_v(array, amplitude=1, frequency=1, shift=0, mask=['ALL']):
    '''Return the sin values of the array.

    Input:
    array - the values in this array will be transform into 2*pi*xi/max_x
            before calculation.

    Output:
    array - the original array
    out_y - an array of sin values of the original array.
    '''
    mask = maskCheck(array, mask)
    array_mod = modx(array)*frequency
    out_y = amplitude*np.sin(array_mod+shift)
    dummy, out_y = maskCalculate(out_y, array, mask)
    return array, out_y


def cos_v(array, amplitude=1, frequency=1, shift=0, mask=['ALL']):
    '''Return the cos values of the array.

    Input:
    array - the values in this array will be transform into 2*pi*xi/max_x
            before calculation.

    Output:
    array - the original array
    out_y - an array of cos values of the original array.
    '''
    mask = maskCheck(array, mask)
    array_mod = modx(array)*frequency
    out_y = amplitude*np.cos(array_mod+shift)
    dummy, out_y = maskCalculate(out_y, array, mask)
    return array, out_y


def sin_sq(array, amplitude, frequency, shift, mask=['ALL']):
    '''Return the sin^2 values of the array.

    Input:
    array - the values in this array will be transform into 2*pi*xi/max_x
            before calculation.

    Output:
    array - the original array
    out_y - an array of sin^2 values of the original array.
    '''
    mask = maskCheck(array, mask)
    array_mod = modx(array)*frequency
    out_y = amplitude*np.square(np.sin(array_mod+shift))
    dummy, out_y = maskCalculate(out_y, array, mask)
    return array, out_y


def cos_sq(array, amplitude, frequency, shift, mask=['ALL']):
    '''Return the cos^2 values of the array.

    Input:
    array - the values in this array will be transform into 2*pi*xi/max_x
            before calculation.

    Output:
    array - the original array
    out_y - an array of cos^2 values of the original array.
    '''
    mask = maskCheck(array, mask)
    array_mod = modx(array)*frequency
    out_y = amplitude*np.square(np.cos(array_mod+shift))
    dummy, out_y = maskCalculate(out_y, array, mask)
    return array, out_y


def offset_v(s_v, start, minimum, direction, fun=None):
    '''Calculate the offsets from a centerline.

    Input:
    s_v - array. Values feed to function to calculate the offsets.
            Can be thought as the x values of the line.
    start - array. Initial values that will be added to offsets.
            Can be thought as the y valuse of the line.
    minimum - minimum value of the offset.
    direction - string. 'left' or 'right'.
                'left' increases; 'right' decreases.
    fun - function used to calculate offsets.

    Output: an array contains offset values.
    '''
    variate = 0
    protect = 0
    if fun is not None:
        dummy, variate = fun(s_v)
        protect = np.amin(variate)
    if direction == "left":
        return variate+minimum+start-protect
    elif direction == "right":
        return -variate-minimum+start+protect


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


def curvature_v(x_v, y_v):
    '''Return an array of curvature.
    
    Inputs:
    x_v - array. x values.
    y_v - array. y values.

    Output:
    s_v - array. curvatures.
    '''
    x1_v = x_v[:-1]

    s_v = slope_v(x_v, y_v)
    curvature = slope_v(x_v, s_v)
    return curvature


def shields(d50, css, s, g_s=1922, g_w=1000):
    ''' Shields equation:
        the default sediment is sand
        data is from https://www.simetric.co.uk/si_materials.htm

        The result has unit in meters.
    '''
    if s == 0 or g_w == 0:
        return -1
    else:
        return (g_s-g_w) * d50 * css / (g_w*s)


def sn_to_xy(x1, y1, x2, y2, n):
    '''Transform a point from sn coodinate system to xy coodinate system.

    Inputs:
    x1, y1 - x, y values of point at the centerline.
    x2, y2 - x, y values of next point at the centerline.
    n - distance from the point we want to transform to the first point at
        centerline.

    Outputs:
    x, y - x, y values of point we want to transform.
    '''
    dx = x2 - x1
    dy = y2 - y1
    ds = math.sqrt(dx**2 + dy**2)
    x = dy/ds*(-n) + x1
    y = dx/ds*n + y1
    return x, y


def sn_to_xy_v(x_v, y_v, n_v):
    '''Transform an array of points from sn coordinate system to xy coordinate system.
    The last point has the same x, y as the second last point to maintain the length.

    Inputs:
    x_v, y_v - x, y values of points at the centerline.
    n_v - offsets from the points to their corresponding centerline points.

    Output:
    out_x, out_y - x, y values of points need to be transformed.
    '''
    x1_v = x_v[:-1]
    x2_v = x_v[1:]
    dx = x2_v - x1_v

    y1_v = y_v[:-1]
    y2_v = y_v[1:]
    dy = y2_v - y1_v

    ds = np.sqrt(np.square(dy) + np.square(dx))
    out_x = dy/ds*(-1*n_v[:-1]) + x1_v
    out_y = dx/ds*(n_v[:-1]) + y1_v
    out_x = np.append(out_x, out_x[-1])
    out_y = np.append(out_y, out_y[-1])
    return out_x, out_y


def xy_to_sn(x1, y1, x2, y2, direction):
    '''Transform a point from xy coodinate system to sn coodinate system.

    Inputs:
    x1, y1 - x, y values of centerline point.
    x2, y2 - x, y values of offset point.
    direction - 'left' indicates positive; 'right' indicates negative.

    Output:
    n - distance from offset point to centerline point.
    '''
    dx = x2 - x1
    dy = y2 - y1
    n = math.sqrt(dx**2 + dy**2)
    if direction == "left":
        return n
    else:
        return -n


def fx(x1, y1, x2, y2, x3):
    '''Linear function calculator.
    
    Inputs:
    x1, y1 - info of point 1.
    x2, y2 - info of point 2.
    x3 - x of point 3.

    Output:
    y of point 3.
    '''
    s = slope(x1, y1, x2, y2)
    dx = x3 - x1
    return s*dx + x1


def pointDist(p1, p2):
    '''Distance from point 1 to point 2.'''
    return math.sqrt((p2[1]-p1[1])**2 + (p2[0]-p1[0])**2)


def dist_v(x_v, y_v):
    '''Distance between neighbors in a line.

        The length will decrease by 1.
    '''
    x1 = x_v[:-1]
    x2 = x_v[1:]
    y1 = y_v[:-1]
    y2 = y_v[1:]

    return (y2-y1)**2 + (x2-x1)**2

def normalDist(slope, point1, point2):
    '''Return distance from point2 to perpendicular line of point1.
    The perpendicular line is in form: as a*x + b*y + c = 0

    slope -- slope at point1
    point1, point2 -- (x, y)
    '''
    if slope == 0:
        return abs(point2[0] - point1[0])

    slope = -1/slope
    a = slope
    b = -1
    c = point1[1] - slope*point1[0]

    return abs(a*point2[0] + b*point2[1] + c) / math.sqrt(a*a + b*b)

def vectorIntersect(v1_p1, v1_p2, v2_p1, v2_p2):
    '''Return boolean if these two vectors intersect with each other

    vi_pj -- j point of vector i
    
    Reference:
    https://stackoverflow.com/questions/563198
    '''
    p = np.array(v1_p1)
    r = np.array(v1_p2) - p
    q = np.array(v2_p1)
    s = np.array(v2_p2) - q

    d = np.cross(r, s)
    if d == 0:
        return False

    t = np.cross((q-p), s) / d
    u = np.cross((q-p), r) / d

    if u <=1 and u >=0:
        return t
    else:
        return False


def linefx(p1, slope, x):
    '''Return y value of x on a given line p1 to p2.'''
    dx = x - p1[0]
    y = p1[1] + dx*slope
    return y


def interpolation(y1, y2, x):
    '''Interpolation function to connect y1, y2 smoothly.
    
    x should be in [0, 1]
    '''
    x = x % 1
    x = 6*x**5 - 15*x**4 + 10*x**3
    return (1-x)*y1 + x*y2


def perlin(x_v, amplitude, wavelength, octave, mask):
    '''Return the result of a perlin function.
    
    x_v - uniformly places vector to calculate perlin function
    amplitude - maximum value of the result
    wavelength - how frequently will we comput a random value,
                 smaller wavelength means higher frequency.
    '''
    mask = maskCheck(x_v, mask)

    if wavelength >= len(x_v):
        wavelength = wavelength % len(x_v)

    final_out = np.array([0]*len(x_v))
    for num in range(octave):
        out = []

        wavelength = int(wavelength)
        y0 = (random()-0.5) * 2
        x0 = 0

        for i in range(wavelength, len(x_v)+1, wavelength):
            x1 = i
            y1 = (random() - 0.5)*2 + y0

            for t in range(x0, x1):
                alpha = (t-x0)/(x1-x0)
                y = interpolation(y0, y1, alpha)
                out.append(y)

            x0 = x1
            y0 = y1

        residue = len(x_v) % wavelength 
        if residue != 0:
            y0 = out[-1]
            y1 = (random() - 0.5)*2 + y0
            for i in range(residue):
                alpha = i / residue
                y = interpolation(y0, y1, alpha)
                out.append(y)

        out = np.array(out)
        out = out*amplitude
        final_out = final_out + out

        amplitude = amplitude/2
        wavelength = wavelength/2

    dummy, final_out = maskCalculate(final_out, x_v, mask)
    return x_v, final_out


def cnoidal(x_v, a, f, s, m, mask):
    '''Return the result of a cnoidal waves.
    
    Inputs:
    x_v - uniformly places array to calculate cnoidal function
    a - amplitude, maximum value of the result
    f - frequency
    s - phase shift
    m - elliptic parameter ([0, 1], closer to 1, sharper peak)
    mask - on and off indication of the function.

    Output:
    x_v, y_v - x, y values of cnoidal wave.
    '''
    mask = maskCheck(x_v, mask)
    dummy, siny = sin_v(x_v, 1, f, s, ['ALL'])
    km = 1/np.sqrt(1-m*siny**2)
    dummy, cosy = cos_v(x_v, 1, f, s, ['ALL'])
    y_v = a*(cosy / 2*km)**2
    ymax = np.amax(y_v)
    y_v = ymax-y_v
    dummy, y_v = maskCalculate(y_v, x_v, mask)
    return x_v, y_v


def step(x_v, a, f, s, m, mask):
    '''Return the result of a wave looked like steps.
    
    x_v - uniformly places vector to calculate perlin function
    a - amplitude, maximum value of the result
    f - frequency
    s - phase shift
    m - elliptic parameter ([0, 1], closer to 1, sharper peak)
    '''
    mask = maskCheck(x_v, mask)
    dummy, siny = sin_v(x_v, 1, f, s, ['ALL'])
    km = 1/np.sqrt(1-m*siny**2)
    dummy, cosy = cos_v(x_v, 1, f, s, ['ALL'])
    out_y = a*(cosy / 2*km)
    dummy, out_y = maskCalculate(out_y, x_v, mask)
    return x_v, out_y


def highCurv(x_v, a, f, ps, p, mask):
    '''
    
    x_v - uniformly places vector to calculate high curvature function.
    a - amplitude, maximum value of the result
    f - frequency
    ps - phase shift
    p - sinuosity (with p>2, we should have high curvature function)
    '''
    m = p*len(x_v)
    w = 2.2*math.sqrt((p-1)/p)
    s = np.array(list(range(int(m))))
    dummy, cosy= cos_v(s, 1, f, ps, '1')
    theta_v = w*cosy
    
    x_pre = 0
    y_pre = 0
    out_x = [x_pre]
    out_y = [y_pre]

    for i in range(len(s)):
        theta = theta_v[i]
        x = x_pre + math.cos(theta)
        y = y_pre + math.sin(theta)

        out_x.append(x)
        out_y.append(y)

        x_pre = x
        y_pre = y

    out_x = np.array(out_x)
    mask = maskCheck(out_x, mask)
    out_y = np.array(out_y)
    out_y_max = np.amax(out_y)
    out_y = out_y/out_y_max * a
    out_x, out_y = maskCalculate(out_y, out_x, mask)
    return out_x, out_y


def unit_vector(vector):
    '''return the unit vector of the vector.'''
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    '''Returns the angle in radians between vectors v1 and v2'''
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def indexBound(x, x_v_sort):
    ind1 = ind2 = len(x_v_sort) - 1
    for i in range(len(x_v_sort)):
        if x_v_sort[i] == x:
            ind1 = ind2 = i
            break
        elif x_v_sort[i] > x:
            ind1 = i - 1
            ind2 = i
            if ind1 < 0:
                ind1 = 0
            break

    return ind1, ind2


#def deleteCycles(line):
#    '''
#    line - list of points, each point is (x, y, ...),
#            the order of points is the trend of the line
#    '''
#    lineGrows_positive = []
#    lineGrows_negative = []
#    circles_index = {}
#    lst_pt = (line[0][0], line[0][1], 0)
#
#    for i in range(1, len(line)):
#        crt_pt = (line[i][0], line[i][1], i)
#        if int(crt_pt[0]) < 0:
#            extend = updateLine(lineGrows_negative, crt_pt, i)
#        else:
#            extend = updateLine(lineGrows_positive, crt_pt, i)
#
#        if extend:
#            lst_pt = crt_pt
#            continue
#
#        ind1, ind2 = detectCycle(lineGrows_positive, lineGrows_negative, line, crt_pt, lst_pt)
#
#        if ind1 != -1:
#            #print(ind1+1, ind2-1)
#            circles_index[ind1+1] = ind2-1
#
#        lst_pt = crt_pt
#
#    if circles_index == {}:
#        return line
#
#    clean_line = []
#    outCircleMarker = [True for i in range(len(line))]
#    for key in circles_index.keys():
#        left = key
#        right = circles_index[key]
#        for i in range(left, right+1):
#            outCircleMarker[i] = False
#
#    clean_line = [line[i] for i in range(len(outCircleMarker)) if outCircleMarker[i]]
#
#    return clean_line


def deleteCycles(line):
    '''
    line - list of points, each point is (x, y, ...),
            the order of points is the trend of the line
    '''
    clean_line = []
    line_copy = [x for x in line]
    pointsDict = {}
    cycleIndex = []

    line_copy = completeLine(line_copy)
    for i in range(0, len(line_copy)):
        crt_pt = line_copy[i]
        crt_pt_xy = (crt_pt[0], crt_pt[1])
        if crt_pt_xy in pointsDict:
            ind1 = pointsDict[crt_pt_xy]
            cycleIndex.append((ind1, i))
        else:
            pointsDict[crt_pt_xy] = i

    line_marker = [True for i in range(len(line_copy))]
    for (left, right) in cycleIndex:
        for i in range(left+1, right):
            line_marker[i] = False

    clean_line = [line_copy[i] for i in range(len(line_copy)) if line_marker[i]]

    return clean_line


def completeLine(ln):
    '''
    1. Round points in ln to int (only x and y values)
    2. Remove duplicate points in ln (if both x and y matches)
    '''
    out_ln = [(round(ln[0][0]), round(ln[0][1]), ln[0][2])]
    lst_pt = out_ln[0]
    for i in range(1, len(ln)):
        crt_pt = (round(ln[i][0]), round(ln[i][1]), ln[i][2])
        lst_pt_xy = (lst_pt[0], lst_pt[1])
        crt_pt_xy = (crt_pt[0], crt_pt[1])

        if lst_pt_xy == crt_pt_xy:
            continue
        
        diff_x = abs(lst_pt_xy[0] - crt_pt_xy[0])
        diff_y = abs(lst_pt_xy[1] - crt_pt_xy[1])

        if diff_x >= 2 or diff_y >= 2:
            steps = int(max(diff_x, diff_y))
            for t in range(1, steps+1):
                alpha = t/steps
                x = round(lst_pt[0]*(1-alpha) + crt_pt[0]*alpha)
                y = round(lst_pt[1]*(1-alpha) + crt_pt[1]*alpha)
                z = lst_pt[2]*(1-alpha) + crt_pt[2]*alpha
                out_ln.append((x, y, z))
        elif diff_x == diff_y:
            out_ln.append((lst_pt[0], crt_pt[1], crt_pt[2]))
            out_ln.append(crt_pt)
        else:
            out_ln.append(crt_pt)

        lst_pt = crt_pt

    return out_ln


def updateLine(linelist, pt, index):
    '''
    linelist - [[(x, y, ind)], ...] x is index of the list
    pt - (x, y)
    '''
    extend = False
    round_x = abs(int(pt[0]))
    if round_x >= len(linelist):
        linelist += [[] for i in range(round_x - len(linelist)+1)] 
        extend = True

    linelist[round_x].append((pt[0], pt[1], index))
    linelist[round_x].sort()

    return extend


#def detectCycle(line_positive, line_negative, orig_line, crt_pt, lst_pt):
#    '''
#    crt_pt - (x, y, index)
#    lst_pt - (x, y, index)
#    '''
# #vectorIntersect(v1_p1, v1_p2, v2_p1, v2_p2):
#    crt_x = int(crt_pt[0])
#    if crt_x < 0 :
#        used_line = line_negative
#    else:
#        used_line = line_positive
#
#    crt_x = abs(crt_x)
#    ind_crt_pt = used_line[crt_x].index(crt_pt)
#    cls_pt = lst_pt
#    ind_cls_pt = ind_crt_pt
#    cls_pt_orig_ind = crt_pt[2] + 1
#    check_ind = cls_pt_orig_ind
#    check_x = crt_x
#
#    while cls_pt_orig_ind == check_ind:
#
#        if used_line[check_x] == []:
#            check_x += 1
#            ind_cls_pt = 0
#            continue
#
#        if ind_cls_pt != len(used_line[check_x]):
#            cls_pt = used_line[check_x][ind_cls_pt]
#            ind_cls_pt += 1
#            cls_pt_orig_ind = cls_pt[2]
#            check_ind -= 1
#        elif check_x != len(used_line) - 1:
#            check_x += 1
#            ind_cls_pt = 0
#        else:
#            return -1, -1
#
#    left_bound = cls_pt_orig_ind - 2
#    if left_bound < 0:
#        left_bound = 0
#    right_bound = cls_pt_orig_ind + 2
#
#    v1_p1 = (lst_pt[0], lst_pt[1])
#    v1_p2 = (crt_pt[0], crt_pt[1])
#    for i in range(left_bound+1, right_bound):
#        v2_p1 = (orig_line[i-1][0], orig_line[i-1][1])
#        v2_p2 = (orig_line[i][0], orig_line[i][1])
#
##        print(v1_p1, v1_p2, v2_p1, v2_p2)
#        if vectorIntersect(v1_p1, v1_p2, v2_p1, v2_p2):
##            print('bound', left_bound+1, right_bound, crt_pt[2])
#            return i-1, crt_pt[2]
#
#    return -1, -1

def perlin2D(x,y,seed=1):
    # permutation table
    #np.random.seed(seed)
    p = np.arange(256,dtype=int)
    np.random.shuffle(p)
    p = np.stack([p,p]).flatten()
    # coordinates of the top-left
    xi = x.astype(int)
    yi = y.astype(int)
    # internal coordinates
    xf = x - xi
    yf = y - yi
    # fade factors
    u = fade(xf)
    v = fade(yf)
    # noise components
    n00 = gradient(p[p[xi]+yi],xf,yf)
    n01 = gradient(p[p[xi]+yi+1],xf,yf-1)
    n11 = gradient(p[p[xi+1]+yi+1],xf-1,yf-1)
    n10 = gradient(p[p[xi+1]+yi],xf-1,yf)
    # combine noises
    x1 = lerp(n00,n10,u)
    x2 = lerp(n01,n11,u) # FIX1: I was using n10 instead of n01
    return lerp(x1,x2,v) # FIX2: I also had to reverse x1 and x2 here

def lerp(a,b,x):
    "linear interpolation"
    return a + x * (b-a)

def fade(t):
    "6t^5 - 15t^4 + 10t^3"
    return 6 * t**5 - 15 * t**4 + 10 * t**3

def gradient(h,x,y):
    "grad converts h to the right gradient vector and return the dot product with (x,y)"
    vectors = np.array([[0,1],[0,-1],[1,0],[-1,0]])
    g = vectors[h%4]
    return g[:,:,0] * x + g[:,:,1] * y
