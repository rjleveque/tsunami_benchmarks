
from pylab import *
from scipy.optimize import least_squares
import matplotlib.animation as animation


def make_corners_fcn(L,phi):
    
    def corners(z):
        """
        convert z = [x0,y0,theta] into a list of corners, assuming (x0,y0) is
        starting corner and moving at angle theta (up from x-axis) to 2nd corner.
        Length of first side is L0, then turn through angle phi1, etc.
        """
        x0,y0,theta = z # unpack
        xc = zeros(len(L)+1)
        yc = zeros(len(L)+1)
        xc[0] = x0
        yc[0] = y0
        phitot = theta
        for k in range(len(L)):
            phitot = phitot + phi[k]
            xc[k+1] = xc[k] + L[k]*cos(phitot)
            yc[k+1] = yc[k] + L[k]*sin(phitot)
        return xc,yc
        
    return corners
            

L0 = 1.
L1 = 1.
L2 = 1.
phi1 = pi/2
phi2 = pi/2
x0 = 0; y0 = 0;
x1 = 0; y1 = L0;
x2 = L1; y2 = L0;
x3 = L1; y3 = 0;
xc_init = array([x0,x1,x2,x3,x0])
yc_init = array([y0,y1,y2,y3,y0])
theta_init = 0.
z_init = array([x0,y0,theta_init])
            
def corners1(z):
    x0,y0,theta = z # unpack
    x1 = x0 + L0*sin(theta)
    y1 = y0 + L0*cos(theta)
    x2 = x1 + L1*sin(theta + phi1)
    y2 = y1 + L1*cos(theta + phi1)
    x3 = x2 + L2*sin(theta + phi1 + phi2)
    y3 = y2 + L2*cos(theta + phi1 + phi2)
    xc = array([x0,x1,x2,x3,x0])
    yc = array([y0,y1,y2,y3,y0])
    return xc,yc
    
L3 = 1
phi3 = pi/2
L = array([L0,L1,L2,L3])
phi = array([0,phi1,phi2,phi3])    
corners = make_corners_fcn(L,phi)
    
def F(z, *args, **kwargs):
    
    # target corners xhat,yhat and corners fcn passed in as args:
    xhat = args[0]
    yhat = args[1]
    corners = args[2]
    
    # compute current guess at corners from z
    xc,yc = corners(z)

    # compute residuals between current corners and target:
    f = zeros(2*len(xc))
    for i in range(len(xc)):
        f[2*i] = xc[i] - xhat[i]
        f[2*i+1] = yc[i] - yhat[i]
    return f
    
    
def remap(xc_hat, yc_hat, z_guess, corners):
    """
    Adjust xc_hat, yc_hat to xc,yc so that original shape is preserved.
    z_guess is initial guess for z defining xc,yc.
    """
    
    result = least_squares(F,z_guess,args=(xc_hat,yc_hat,corners))
    xc,yc = corners(result['x'])
    theta = result['x'][2]
    return xc,yc,theta
    
def velocities_shear():
    u = lambda x,y,t: y
    v = lambda x,y,t: zeros(x.shape)
    return u,v
    
def make_corner_paths(L,phi,z0,u,v,t0,dt,nsteps):
    """
    make a list of lists, one for each time t0, t0+dt, ... t0 + nsteps*dt
    """
    corners = make_corners_fcn(L,phi)
    xc,yc = corners(z0)
    ncorners = len(xc)
    uc = [u(xc[k],yc[k],t0) for k in range(ncorners)]
    vc = [v(xc[k],yc[k],t0) for k in range(ncorners)]
    corner_paths = [[t0, vstack([xc,yc,uc,vc]).T]]
    z_np1 = z0
    for n in range(nsteps):
        xc_hat = []
        yc_hat = []
        t_n,corners_n = corner_paths[-1] # from previous time step
        t_np1 = t_n + dt
        for k in range(ncorners):
            xc_n, yc_n, uc_n, vc_n = corners_n[k,:]  # unpack k'th row
            xc_hat.append(xc_n + dt*uc_n)
            yc_hat.append(yc_n + dt*vc_n)
        # remap to original shape:
        xc_np1, yc_np1, theta_np1 = remap(xc_hat, yc_hat, z_np1, corners)
        z_np1 = (xc_np1[0], yc_np1[0], theta_np1)
        uc_np1 = u(xc_np1, yc_np1, t_np1)
        vc_np1 = v(xc_np1, yc_np1, t_np1)
        corners_np1 = vstack([xc_np1,yc_np1,uc_np1,vc_np1]).T
        corner_paths.append([t_np1, corners_np1])
    
    return corner_paths      
            
        
def test_corner_paths_square():
    u,v = velocities_shear()
    t0 = 0.
    nsteps = 8
    dt = 2*pi/nsteps
    z0 = [3,1,pi/4]
    corner_paths = make_corner_paths(L,phi,z0,u,v,t0,dt,nsteps)
    
    figure(1);clf();
    for k in range(8):
        tk,cpk = corner_paths[k]
        plot(cpk[:,0],cpk[:,1],'b')
        # one edge different color to show orientation:
        plot(cpk[:2,0],cpk[:2,1],'c')
    axis('square')
    
    return corner_paths

def test_corner_paths_triangle():
    L = [1,1,1]
    phi = [2*pi/3, 2*pi/3, 2*pi/3]
    z0 = [0,0.5,0]
    
    u,v = velocities_shear()
    t0 = 0.
    nsteps = 8
    dt = 2*pi/nsteps
    
    corner_paths = make_corner_paths(L,phi,z0,u,v,t0,dt,nsteps)
    
    figure(1);clf();
    for k in range(nsteps):
        tk,cpk = corner_paths[k]
        plot(cpk[:,0],cpk[:,1],'b')
        # one edge different color to show orientation:
        plot(cpk[:2,0],cpk[:2,1],'c')
    axis('square')
    
    return corner_paths
