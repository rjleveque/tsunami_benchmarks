"""
Make piecewise linear topography for wave tank.
"""

from pylab import *
from clawpack.geoclaw_1d import nonuniform_grid_tools

xlower = 0.75
xupper = 43.75
B0 = -0.9056

xzpairs = [(xlower, B0),       # left edge
           (10, B0),           # start of first slope
           (17.5, 0.5+B0),     # start of beach
           (32.5, 1+B0),       # on shore
           (xupper, 1+B0)]     # right edge

topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

mx = 1000
hmin = 0.05  # use uniform grid in shallower water
#hmin = 5  # try forcing a uniform grid everywhere
#hmin = 1  # try uniform grid near start of shelf

nonuniform_grid_tools.make_celledges_cfl(xlower, xupper, mx, topo_fcn,
        hmin, fname='celledges.data', plot_topo=True)
