
from clawpack.geoclaw import topotools
from pylab import *

dx=.0025;

def maketopo():
    """
    Output topography files
    """
    nxpoints = 201
    nypoints = 21
    xlower = 0.
    xupper = 20.
    ylower = 0.
    yupper = 1.52
    outfile= "domain.tt1"
    topotools.topo1writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

    nxpoints = 201
    nypoints = 201
    xlower = 4.5
    xupper = 5.5
    ylower = 0.25
    yupper = 1.25
    outfile= "hump.tt1"
    topotools.topo1writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)


def topo(x,y):
    """
    flat with hump
    """

    z = zeros(x.shape)

    # island addition
    # from paper
    xo=5.0;
    yo=1.52/2;
    slope_d=8;
    base_r=0.75/2;
    top_r=0.05/2;
    hi=0.049;

    dist = sqrt( (x-xo)**2+(y-yo)**2 );
    zz = hi*(1-(dist-top_r)/(base_r-top_r))
    zz = where(zz<hi, zz, hi)
    z = where(dist<base_r, zz, 0.)

    z = z - 0.054
    
    return z

def plot_topo():

    topo1 = topotools.Topography()
    topo1.read('domain.tt1',1)
    figure(figsize=(12,4))
    ax = axes()
    topo1.plot(axes=ax,add_colorbar=False)

    topo2 = topotools.Topography()
    topo2.read('hump.tt1',1)
    topo2.plot(axes=ax)
    ax.set_xlim(0,9.84)
    ax.set_ylim(0,1.52)
    axis('scaled')

if __name__=='__main__':
    maketopo()
    
