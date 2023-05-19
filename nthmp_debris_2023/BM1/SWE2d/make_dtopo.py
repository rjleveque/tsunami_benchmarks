
from pylab import *
from clawpack.geoclaw import dtopotools

dtopo = dtopotools.DTopography()


rad = 0.3
xcentroid = 34.94
ycentroid = 0.82

x1b = xcentroid - 2*rad
x2b = xcentroid + 2*rad
y1b = ycentroid - 2*rad
y2b = ycentroid + 2*rad

x = linspace(x1b-1*rad,x2b+1*rad,301)
y = linspace(y1b-1*rad,y2b+1*rad,301)

X,Y = meshgrid(x,y)
dZb = where(logical_and(logical_and(X>x1b,X<x2b), logical_and(Y>y1b,Y<y2b)),
            1., 0.)
            #0.051, 0.)

dZ0 = zeros(X.shape)
times = array([33,34,35,35.3])
dZ = empty((len(times), len(y), len(x)))
dZ[0,:,:] = dZ0
dZ[1,:,:] = dZb
dZ[2,:,:] = dZb
dZ[3,:,:] = dZ0

dtopo.x = x
dtopo.y = y
dtopo.X = X
dtopo.Y = Y
dtopo.dZ = dZ
dtopo.times = times

fname = 'dtopo4.dtt1'
dtopo.write(fname, dtopo_type=1)
print('Created ',fname)

#============================

xcentroid = 34.64
ycentroid = 0.52

x1b = xcentroid - 1*rad
x2b = xcentroid + 1*rad
y1b = ycentroid - 1*rad
y2b = ycentroid + 1*rad

x = linspace(x1b-1*rad,x2b+1*rad,201)
y = linspace(y1b-1*rad,y2b+1*rad,201)

X,Y = meshgrid(x,y)
dZb = where(logical_and(logical_and(X>x1b,X<x2b), logical_and(Y>y1b,Y<y2b)),
            1., 0.)
            #0.051, 0.)

dZ0 = zeros(X.shape)
times = array([33,34,35,35.3])
dZ = empty((len(times), len(y), len(x)))
dZ[0,:,:] = dZ0
dZ[1,:,:] = dZb
dZ[2,:,:] = dZb
dZ[3,:,:] = dZ0

dtopo.x = x
dtopo.y = y
dtopo.X = X
dtopo.Y = Y
dtopo.dZ = dZ
dtopo.times = times

fname = 'dtopo1.dtt1'
dtopo.write(fname, dtopo_type=1)
print('Created ',fname)
