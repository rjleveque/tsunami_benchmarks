
from pylab import *
from clawpack.geoclaw import topotools

wavetank = topotools.Topography()
wavetank.read('wavetank.tt2',2)
blocks = topotools.Topography()
blocks.read('blocks.tt2',2)
fig,ax = subplots(figsize=(12,3))

wavetank.plot(axes=ax,limits=(-1,1))
blocks.plot(axes=ax,limits=(-1,1),add_colorbar=False)
axis('equal')
axis([0,41.29,-4,4]) 

