
from pylab import *
from clawpack.geoclaw import dtopotools

dtopo = dtopotools.DTopography()

# centers of circles in rectangular grid:
dd = 0.153 # distance between centers
#xg1 = arange(31.29+0.5*dd, 31.29+4*dd, dd)
#yg1 = arange(-2*dd, 2.1*dd, dd)

x1b = 31.29
x2b = x1b + 4*dd
y1b = -2*dd
y2b = 2*dd

x = linspace(x1b-2*dd,x2b+2*dd,401)
y = linspace(y1b-1*dd,y2b+1*dd,301)

X,Y = meshgrid(x,y)
dZb = where(logical_and(logical_and(X>x1b,X<x2b), logical_and(Y>y1b,Y<y2b)),
            1., 0.)
            #0.051, 0.)

dZ0 = zeros(X.shape)
times = array([27,28,29,30,31,32])
dZ = empty((len(times), len(y), len(x)))
dZ[0,:,:] = dZ0
dZ[1,:,:] = dZb
dZ[2,:,:] = dZb
dZ[3,:,:] = dZb
dZ[4,:,:] = dZb
dZ[5,:,:] = dZ0

dtopo.x = x
dtopo.y = y
dtopo.X = X
dtopo.Y = Y
dtopo.dZ = dZ
dtopo.times = times

fname = 'dtopo.dtt1'
dtopo.write(fname, dtopo_type=1)
print('Created ',fname)
