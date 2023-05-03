"""
Make piecewise linear topography for wave tank.
"""

from pylab import *
from clawpack.geoclaw_1d import nonuniform_grid_tools

xlower = 0.
xupper = 41.29
B0 = -0.87

xzpairs = [(xlower, B0),       # left edge
           (11.29, B0),        # start of slope
           (31.29, 1+B0),      # end of slope
           (xupper, 1+B0)]     # right edge

topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

mx = 1000
hmin = 0.05  # use uniform grid in shallower water
#hmin = 5  # try forcing a uniform grid everywhere
#hmin = 1  # try uniform grid near start of shelf

nonuniform_grid_tools.make_celledges_cfl(xlower, xupper, mx, topo_fcn,
        hmin, fname='celledges.data', plot_topo=True)
