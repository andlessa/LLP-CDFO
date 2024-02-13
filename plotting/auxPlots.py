
import numpy as np
import itertools
from scipy.interpolate import  griddata
from matplotlib import pyplot as plt


def interpolateData(x,y,z,nx=200,ny=200,method='linear',fill_value=np.nan,xnew=None,ynew=None):

    if x.min() == x.max() or y.min() == y.max(): # Can not interpolate
        return None,None,None
    elif xnew is None or ynew is None:
        xnew = np.linspace(x.min(),x.max(),nx)
        ynew = np.linspace(y.min(),y.max(),ny)

    xi = np.array([list(v) for v in itertools.product(xnew,ynew)])
    znew = griddata(list(zip(x,y)),z,xi=xi, 
                    method=method,fill_value=fill_value)
    znew = np.reshape(znew,(len(xnew),len(ynew)))
    xnew,ynew  = np.meshgrid(xnew,ynew,indexing='ij')

    return xnew,ynew,znew

def getContours(x,y,z,contourValues):
    

    contours = plt.contour(x, y, z, contourValues)
    plt.close()

    contoursDict = {}

    for i,item in enumerate(contours.collections):
        cV = contourValues[i]
        xData = []
        yData = []
        for p in item.get_paths():
            v = p.vertices
            xData += list(v[:, 0])
            yData += list(v[:, 1])
        if len(xData) == 0:
            continue
        contoursDict[cV] = np.array(list(zip(xData,yData)))
    
    return contoursDict

def saveContours(contoursDict,fname,header):

    with open(fname,'w') as f:
        for cV,data in contoursDict.items():
            np.savetxt(f,data,fmt='%.4e',delimiter=',',header=header,comments='\n\n# Contour value=%1.2f \n' %cV)
    print('Contours saved to %s' %fname)

def readContours(fname):

    contoursDict = {}
    with open(fname,'r') as f:
        dataBlocks = f.read().split('#')[1:]
        for data in dataBlocks: 
            data = data.splitlines()
            cV = eval(data[0].split('=')[1])
            dataPts = np.genfromtxt(data,delimiter=',',names=True,skip_header=1)
            contoursDict[cV] = dataPts

    return contoursDict

def label_line(fig,line, label_text, 
               near_i=None, near_x=None, near_y=None, 
               rotation_offset=0, offset=(0,0),fontsize=13,
               xmin=None,rotation=None,boxalpha=1.0):
    """call 
        l, = plt.loglog(x, y)
        label_line(l, "text", near_x=0.32)
    """
    def put_label(i):
        """put label at given index"""
        i = min(i, len(x)-2)
        dx = sx[i+1] - sx[i]
        dy = sy[i+1] - sy[i]
        if rotation is None:
            rot = np.rad2deg(np.arctan2(dy, dx)) + rotation_offset
        else:
            rot = rotation
        pos = [(x[i] + x[i+1])/2. + offset[0], (y[i] + y[i+1])/2 + offset[1]]
        if pos[0] > xmin:
            plt.text(pos[0], pos[1], label_text, size=fontsize, 
                     rotation=rot, color = line.get_color(),
                     ha="center", va="center", bbox = dict(ec='1',fc='1',alpha=boxalpha))

    x = line.get_xdata()
    y = line.get_ydata()
    ax = fig.get_axes()[0]
    if ax.get_xscale() == 'log':
        sx = np.log10(x)    # screen space
    else:
        sx = x
    if ax.get_yscale() == 'log':
        sy = np.log10(y)
    else:
        sy = y

    # find index
    if near_i is not None:
        i = near_i
        if i < 0: # sanitize negative i
            i = len(x) + i
        put_label(i)
    elif near_x is not None:
        for i in range(len(x)-2):
            if (x[i] < near_x and x[i+1] >= near_x) or (x[i+1] < near_x and x[i] >= near_x):
                put_label(i)
    elif near_y is not None:
        for i in range(len(y)-2):
            if (y[i] < near_y and y[i+1] >= near_y) or (y[i+1] < near_y and y[i] >= near_y):
                put_label(i)
    else:
        raise ValueError("Need one of near_i, near_x, near_y")