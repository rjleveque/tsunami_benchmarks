"""
Create the BM1_leveque.txt file requested by Pat Lynett.
"""

from pylab import *
from clawpack.visclaw.data import ClawPlotData

plotdata = ClawPlotData()

toffset = 0.

if 1:
    plotdata.outdir = '_output_manning010_cfl090'
    fname = 'BM1_leveque_1.txt'
if 0:
    plotdata.outdir = '_output_manning010_cfl089'
    fname = 'BM1_leveque_2.txt'
if 0:
    plotdata.outdir = '_output_manning015_cfl090'
    fname = 'BM1_leveque_3.txt'
if 0:
    plotdata.outdir = '_output_manning015_cfl089'
    fname = 'BM1_leveque_4.txt'

figure(figsize=(8,12))
clf()

# ---  Gauge 1 ---

d = loadtxt('s1u.txt')
t1u = d[:,0]
s1u = d[:,1]

d = loadtxt('s1v.txt')
t1v = d[:,0]
s1v = d[:,1]

g = plotdata.getgauge(1)
u1 = g.q[1,:] / g.q[0,:]
v1 = g.q[2,:] / g.q[0,:]
t1 = g.t

subplot(4,1,1)
plot(t1u+toffset,s1u,'b',label='Experiment')
plot(t1, u1, 'r',label='GeoClaw')
ylabel('u (m/s)')
legend(loc='upper right')
title(plotdata.outdir)

subplot(4,1,2)
plot(t1v+toffset,s1v,'b',label='Experiment')
plot(t1, v1, 'r',label='GeoClaw')
ylabel('v (m/s)')

# ---  Gauge 2 ---

d = loadtxt('s2u.txt')
t2u = d[:,0]
s2u = d[:,1]

d = loadtxt('s2v.txt')
t2v = d[:,0]
s2v = d[:,1]

g = plotdata.getgauge(2)
u2 = g.q[1,:] / g.q[0,:]
v2 = g.q[2,:] / g.q[0,:]
t2 = g.t



subplot(4,1,3)
plot(t2u+toffset,s2u,'b',label='Experiment')
plot(t2, u2, 'r',label='GeoClaw')
ylabel('u (m/s)')
legend(loc='upper right')

subplot(4,1,4)
plot(t2v+toffset,s2v,'b',label='Experiment')
plot(t2, v2, 'r',label='GeoClaw')
ylabel('v (m/s)')

show()

gdata = vstack([t1, u1, v1, u2, v2]).T
savetxt(fname,  gdata, fmt='%.12e')

print "Created ",fname
