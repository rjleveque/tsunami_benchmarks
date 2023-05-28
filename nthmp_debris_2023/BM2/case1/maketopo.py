
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw import topotools
from pylab import *

def maketopo_blocks():
    x1b1 = 35.29
    x2b1 = x1b1 + 0.4
    y1b1 = -0.6
    y2b1 = y1b1 + 0.4
    
    x1b2 = 35.29
    x2b2 = x1b2 + 0.4
    y1b2 = 0.2
    y2b2 = y1b2 + 0.4
    
    B0 = -0.87
    zb = 0.3  # height of blocks
    
    x = linspace(x1b1-0.1, x2b1+0.1, 101)
    y = linspace(y1b1-0.1, y2b2+0.1, 101)
    X,Y = meshgrid(x,y)
    Z = where(logical_and(logical_and(X>=x1b1, X<=x2b1),
                          logical_and(Y>=y1b1, Y<=y2b1)),  
                          1+B0+zb, 1+B0)
    
    Z = where(logical_and(logical_and(X>=x1b2, X<=x2b2),
                          logical_and(Y>=y1b2, Y<=y2b2)),  
                          1+B0+zb, Z)
    
    topo = topotools.Topography()
    topo.set_xyZ(x,y,Z)
    
    topo.write('blocks.tt2', Z_format='%.3f', topo_type=2)


def maketopo_pwlinear():
    """
    Output topography file for the entire domain
    """
    nxpoints = 501
    nypoints = 5
    xlower = 0.e0
    xupper = 41.29
    ylower = -20.e0
    yupper = 20.e0
    outfile= "wavetank.tt2"     
    topotools.topo2writer(outfile,topo_pwlinear,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def topo_pwlinear(x,y):
    """
    piecewise linear
    """
    
    z = zeros(x.shape)
    z = where(x<11.29, z, (x-11.29)/20.)
    z = where(x<31.29, z, 1.)
    z = z - 0.87  # adjust so sea level at 0
    return z


if __name__=='__main__':
    #maketopo_blocks()
    maketopo_pwlinear()
