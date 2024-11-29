
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

def make_corner_paths(debris,u,v,t0,dt,nsteps):
    """
    make a list of lists, one for each time t0, t0+dt, ... t0 + nsteps*dt
    """
    corners = make_corners_fcn(L,phi)
    xc,yc = debris.get_corners(debris.z0)
    ncorners = len(xc)
    uc = [u(xc[k],yc[k],t0) for k in range(ncorners)]
    vc = [v(xc[k],yc[k],t0) for k in range(ncorners)]
    corner_paths = [[t0, vstack([xc,yc,uc,vc]).T]]
    z_np1 = debris.z0
    for n in range(nsteps):
        xc_hat = []
        yc_hat = []
        t_n,corners_n = corner_paths[-1] # from previous time step
        t_np1 = t_n + dt
        for k in range(ncorners):
            xk_n, yk_n, uk_n, vk_n = corners_n[k,:]  # unpack k'th row
            xc_hat.append(xk_n + dt*uk_n)
            yc_hat.append(yk_n + dt*vk_n)
        # remap to original shape:
        xc_np1, yc_np1, theta_np1 = remap(xc_hat, yc_hat, z_np1, corners)
        z_np1 = (xc_np1[0], yc_np1[0], theta_np1)
        uc_np1 = u(xc_np1, yc_np1, t_np1)
        vc_np1 = v(xc_np1, yc_np1, t_np1)
        corners_np1 = vstack([xc_np1,yc_np1,uc_np1,vc_np1]).T
        corner_paths.append([t_np1, corners_np1])
    
    return corner_paths      

        
def test_corner_paths_square():
    
    L = [1,1,1,1]
    phi = [pi/2, pi/2, pi/2, pi/2]
    z0 = [0,0.5,0]
    
    u,v = velocities_shear()
    t0 = 0.
    nsteps = 8
    dt = 2*pi/nsteps
    z0 = [3,1,pi/4]
    corner_paths = make_corner_paths(L,phi,z0,u,v,t0,dt,nsteps)
    
    figure(1);clf();
    for k in range(nsteps):
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


# ===============================

class DebrisObject():
    
    def __init__(self):
        
        self.L = [1.]  # list of edge lengths between corners
        self.phi = [0.]  # list of turning angles at each corner
        self.rho = 0.  # density (g/m^3)
        self.bottom_area = 1.  # area in contact with bottom for friction (m^2)
        self.face_width = 1.   # cross-section area hit by flow (m)
        self.height = 1.  # height of object
        self.advect = True  # passively advected?
        self.friction_static = 0.2
        self.friction_kinetic = 0.1
        self.rho_water = 1000. # density of water (g/m^3)
        self._mass = None # computed from other parameters
        
        self.z0 = [0.,0.,0.]  # location and orientation [x0,y0,theta0]
        
        # these can be generated using velocity field from CFD simulation:
        self.times = []
        self.corner_paths = [] # should be list of [t_n, corners_n] where
                               # corners_n is array with columns [xc,yc,uc,vc]
                               # at each time t_n in self.times
    @property
    def mass(self):
        if self._mass is None:
            self._mass = self.rho * self.bottom_area * self.height
        return self._mass
        
    def get_corners(self, z):    
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

def make_corner_paths_accel(debris,h,u,v,t0,dt,nsteps):
    """
    make a list of lists, one for each time t0, t0+dt, ... t0 + nsteps*dt
    """
    
    # initial corner locations:
    xc,yc = debris.get_corners(debris.z0)
    ncorners = len(xc)
    
    # initial velocities:
    uc = zeros(ncorners)
    vc = zeros(ncorners)
    
    corner_paths = [[t0, vstack([xc,yc,uc,vc]).T]]
    z_np1 = debris.z0
    for n in range(nsteps):
        xc_hat = []
        yc_hat = []
        uc_np1 = []
        vc_np1 = []
        t_n,corners_n = corner_paths[-1] # from previous time step
        t_np1 = t_n + dt
        
        # move corners based on velocities uc_n, vc_n computed at end of last
        # step, but call these uc_hat, vc_hat since they don't maintain rigid
        # body constraint yet:
        for k in range(ncorners):
            xk_n, yk_n, uk_n, vk_n = corners_n[k,:]  # unpack k'th row
            xc_hat.append(xk_n + dt*uk_n)
            yc_hat.append(yk_n + dt*vk_n)

        # remap to original shape, maintaining rigidity:
        xc_np1, yc_np1, theta_np1 = remap(xc_hat, yc_hat, z_np1, corners)
        # corresponding z vector for new position:
        z_np1 = (xc_np1[0], yc_np1[0], theta_np1)
    
        # compute new velocity at n+1 after remapping:
        
        for k in range(ncorners):

            xk_np1 = xc_np1[k]
            yk_np1 = yc_np1[k]
            
            # need uk_n, vk_n to update using accel, values from previous step:
            xk_n, yk_n, uk_n, vk_n = corners_n[k,:]
            # recompute based on actual distance moved, after remapping:
            uk_n = (xk_np1 - xk_n)/dt
            vk_n = (yk_np1 - yk_n)/dt

            # fluid velocities at this corner:
            uk_f = u(xk_np1, yk_np1, t_np1)
            vk_f = v(xk_np1, yk_np1, t_np1)
            
            if debris.advect:
                # to advect with flow, corner velocity = fluid velocity:
                uk_np1 = uk_f
                vk_np1 = vk_f
            else:
                du = uk_f - uk_n
                dv = vk_f - vk_n
                sk_f = sqrt(du**2 + dv**2)
                Ffluid_x = 0.5 * debris.rho_water * sk_f * du * face_area
                Ffluid_y = 0.5 * debris.rho_water * sk_f * dv * face_area
                #Ffluid_x = 0.
                #Ffluid_y = 0.
                uk_np1 = uk_n + dt*Ffluid_x / self.corner_mass
                vk_np1 = vk_n + dt*Ffluid_y / self.corner_mass
                
            uc_np1.append(uk_np1)
            vc_np1.append(vk_np1)
        
        corners_np1 = vstack([xc_np1,yc_np1,uc_np1,vc_np1]).T
        corner_paths.append([t_np1, corners_np1])
    
    return corner_paths         
    
def test_corner_paths_accel():
    L = [1,1,1,1]
    phi = [pi/2, pi/2, pi/2, pi/2]
    z0 = [0,0.5,0]
    
    debris = DebrisObject()
    debris.rho = 500.
        
    u,v = velocities_shear()
    h = lambda x,y,t: 10.  # fluid depth
    
    t0 = 0.
    nsteps = 20
    dt = 2*pi/nsteps
    z0 = [3,1,pi/4]
    corner_paths_a = make_corner_paths_accel(debris,h,u,v,t0,dt,nsteps)
    corner_paths = make_corner_paths(debris,u,v,t0,dt,nsteps)
    
    figure(1);clf();
    for k in range(nsteps):
        tk,cpk = corner_paths[k]
        plot(cpk[:,0],cpk[:,1],'b')
        # one edge different color to show orientation:
        #plot(cpk[:2,0],cpk[:2,1],'c')
        tk,cpk = corner_paths_a[k]
        plot(cpk[:,0],cpk[:,1],'r')
    axis('square')
    
    return corner_paths, corner_paths_a
