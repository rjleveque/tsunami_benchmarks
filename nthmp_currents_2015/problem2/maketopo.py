
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw import topotools
from pylab import *

def maketopo_hilo():

    x = loadtxt('x.txt')
    y = loadtxt('y.txt')
    z = loadtxt('z.txt')
    
    # modify x and y so that cell size is truly uniform:
    dx = 1. / (3.*3600.)   # 1/3"
    xx = linspace(x[0], x[-1], len(x))
    yy = linspace(y[-1], y[0], len(y))
    zz = flipud(z)
    
    topo = topotools.Topography()
    topo.x = xx
    topo.y = yy
    topo.Z = zz
    
    topo.write('hilo_flattened.tt2',topo_type=2)

    

def maketopo_flat():
    """
    Output topography file for the entire domain
    """
    nxpoints = 201
    nypoints = 301
    xlower = 204.812
    xupper = 205.012
    ylower = 19.7
    yupper = 20.0
    outfile= "flat.tt2"     
    topotools.topo2writer(outfile,topo_flat,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def topo_flat(x,y):
    """
    flat
    """
    z = where(x < 204.91213, 30., -30.)
    return z

def plot_topo_big():
    figure(figsize=(8,12))
    topo1 = topotools.Topography()
    topo1.read('flat.tt2',2)
    contourf(topo1.x,topo1.y,topo1.Z,linspace(-30,20,51), extend='both')
    topo2 = topotools.Topography()
    topo2.read('hilo_flattened.tt2',2)
    contourf(topo2.x,topo2.y,topo2.Z,linspace(-30,20,51), extend='both')
    x1 = 204.90028
    x2 = 204.96509
    y1 = 19.71
    y2 = 19.95
    plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],'w')
    axis('scaled')
    colorbar()

def plot_topo():
    figure(figsize=(12,8))
    topo1 = topotools.Topography()
    topo1.read('flat.tt2',2)
    contourf(topo1.x,topo1.y,topo1.Z,linspace(-30,20,51), extend='both')
    topo2 = topotools.Topography()
    topo2.read('hilo_flattened.tt2',2)
    contourf(topo2.x,topo2.y,topo2.Z,linspace(-30,20,51), extend='both')
    colorbar()
    x1 = 204.9
    x2 = 204.955
    y1 = 19.715
    y2 = 19.755
    axis([x1,x2,y1,y2])
    gca().set_aspect(1./cos(y1*pi/180.))
    ticklabel_format(format='plain',useOffset=False)
    contour(topo2.x,topo2.y,topo2.Z,[0.],colors='k')
    plot([204.9447],[19.7308], 'ko')  # from BM description
    plot([204.9437],[19.7307], 'ro')  # closer to pier

    # from <http://tidesandcurrents.noaa.gov/stationhome.html?id=1617760>
    # location is listed as: 19 degrees 43.8' N,  155 degrees, 3.3' W
    xg = 360 - (155 + 3.3/60.)
    yg = 19 + 43.8/60.
    plot([xg],[yg], 'bo')  

    #gauges.append([1125, 204.91802, 19.74517, 0., 1.e9]) #Hilo
    #gauges.append([1126, 204.93003, 19.74167, 0., 1.e9]) #Hilo
    #gauges.append([3333, 204.93, 19.7576,  0., 1.e9]) 


if __name__=='__main__':
    maketopo_hilo()
    maketopo_flat()
