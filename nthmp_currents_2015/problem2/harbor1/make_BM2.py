"""
Create the BM2 files requested by Pat Lynett.
"""

from pylab import *
from scipy import interpolate

from clawpack.visclaw.data import ClawPlotData

plotdata = ClawPlotData()

plotdata.outdir = '_output_1-3sec_alltime'

#tfinal = 4.9 * 3600.
tfinal = 6.4 * 3600.  # for alltime
dt = 1.  # time increment for output files
tout = arange(0., tfinal, dt)   

g = plotdata.getgauge(3333)
p = interpolate.interp1d(g.t, g.q[3,:]) # interpolate surface
g3333_eta = p(tout)

g = plotdata.getgauge(7761)
p = interpolate.interp1d(g.t, g.q[3,:]) # interpolate surface
g7761_eta = p(tout)

g = plotdata.getgauge(1125)
u = g.q[1,:]/g.q[0,:]
v = g.q[2,:]/g.q[0,:]
s = sqrt(u**2 + v**2)
p = interpolate.interp1d(g.t, s) # interpolate speed
g1125_speed = p(tout)

g = plotdata.getgauge(1126)
u = g.q[1,:]/g.q[0,:]
v = g.q[2,:]/g.q[0,:]
s = sqrt(u**2 + v**2)
p = interpolate.interp1d(g.t, s) # interpolate speed
g1126_speed = p(tout)


d = vstack([tout, g3333_eta, g7761_eta, g1125_speed, g1126_speed]).T

savetxt('BM2_leveque.txt',d,fmt='%.12e ')

