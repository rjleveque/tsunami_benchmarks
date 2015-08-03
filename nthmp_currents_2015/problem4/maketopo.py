
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw import topotools
from pylab import *

def maketopo_onshore():

    try:
        x = loadtxt('x_onshore.txt')
        y = loadtxt('y_onshore.txt')
        z = loadtxt('z_onshore.txt')
    except:
        raise Exception("Did you create onshore topo files?  See README.md")
    
    # modify x and y so that cell size is truly uniform:
    xx = arange(x.min(),x.min() + 0.01*(len(x)-1) + .001, 0.01)
    yy = arange(y.min(),y.min() + 0.01*(len(y)-1) + .001, 0.01)
    zz = z - 0.97   # shift so sea_level = 0
    
    topo = topotools.Topography()
    topo.x = xx
    topo.y = yy
    topo.Z = zz
    
    topo.write('seaside_onshore.tt2',topo_type=2)

    figure()
    contourf(xx,yy,zz.T,linspace(-1,0.4,15))
    axis('scaled')
    colorbar()
    

def maketopo_pwlinear():
    """
    Output topography file for the entire domain
    """
    nxpoints = 501
    nypoints = 5
    xlower = 0.e0
    xupper = 50.e0
    ylower = -20.e0
    yupper = 20.e0
    outfile= "pwlinear2.topotype1"     
    topotools.topo1writer(outfile,topo_pwlinear,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def topo_pwlinear(x,y):
    """
    piecewise linear
    """
    
    z = zeros(x.shape)
    z = where(x<10., z, (x-10)/15.)
    z = where(x<17.5, z, 0.5 + (x-17.5)/30.)
    z = where(x<32.5, z, 1.0) 
    z = z - 0.97  # adjust so sea level at 0
    return z


if __name__=='__main__':
    maketopo_onshore()
    maketopo_pwlinear()
