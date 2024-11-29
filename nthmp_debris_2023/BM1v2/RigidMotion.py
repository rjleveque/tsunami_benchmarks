
from pylab import *
from scipy.optimize import least_squares
import matplotlib.animation as animation
from clawpack.visclaw import animation_tools, plottools, geoplot


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
    
    def F(z, *args, **kwargs):
        # Function for least squares fitting
        # compute current guess at corners from z:
        xc,yc = corners(z)

        # compute residuals between current corners and target:
        f = zeros(2*len(xc))
        for i in range(len(xc)):
            f[2*i] = xc[i] - xc_hat[i]
            f[2*i+1] = yc[i] - yc_hat[i]
        return f
    
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
    info = {'friction': None}
    corner_paths = [[t0, vstack([xc,yc,uc,vc]).T, info]]
    z_np1 = debris.z0
    for n in range(nsteps):
        xc_hat = []
        yc_hat = []
        t_n,corners_n,info = corner_paths[-1] # from previous time step
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
        tk,cpk,info = corner_paths[k]
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
        tk,cpk,info = corner_paths[k]
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
        self.grav = 9.81  # gravitational acceleration
        
        self.z0 = [0.,0.,0.]  # location and orientation [x0,y0,theta0]
        
        # these can be generated using velocity field from CFD simulation:
        self.times = []
        self.corner_paths = [] # should be list of [t_n, corners_n, info] where
                               # corners_n is array with columns [xc,yc,uc,vc]
                               # at each time t_n in self.times
                               # and info is a dictionary for info to save

    @property
    def mass(self):
        return self.rho * self.bottom_area * self.height
        
    @property
    def draft(self):
        return (self.rho / self.rho_water) * self.height
                
    def get_corners(self, z):    
        """
        convert z = [x0,y0,theta] into a list of corners, assuming (x0,y0) is
        starting corner and moving at angle theta (up from x-axis) to 2nd corner.
        Length of first side is L0, then turn through angle phi1, etc.
        """
        x0,y0,theta = z # unpack
        xc = zeros(len(self.L)+1)
        yc = zeros(len(self.L)+1)
        xc[0] = x0
        yc[0] = y0
        phitot = theta
        #import pdb; pdb.set_trace()
        for k in range(len(self.L)):
            phitot = phitot + self.phi[k]
            xc[k+1] = xc[k] + self.L[k]*cos(phitot)
            yc[k+1] = yc[k] + self.L[k]*sin(phitot)
        return xc,yc

def make_corner_paths_accel(debris,h,u,v,t0,dt,nsteps,verbose=False):
    """
    make a list of lists, one for each time t0, t0+dt, ... t0 + nsteps*dt
    """
    
    # initial corner locations:
    xc,yc = debris.get_corners(debris.z0)
    ncorners = len(xc)
    
    # initial velocities:
    uc = zeros(ncorners)
    vc = zeros(ncorners)
    
    info = {'friction': 'static'}  # since not moving at t0
    corner_paths = [[t0, vstack([xc,yc,uc,vc]).T, info]]
    z_np1 = debris.z0
    
    for n in range(nsteps):
        xc_hat = []
        yc_hat = []
        #uc_np1 = []
        #vc_np1 = []
        uc_np1 = zeros(ncorners)  # may get modified below
        vc_np1 = zeros(ncorners)
        t_n,corners_n,info_n = corner_paths[-1] # from previous time step
        t_np1 = t_n + dt
        #info_np1 = info_n.copy()  
        
        # move corners based on velocities uc_n, vc_n computed at end of last
        # step, but call these uc_hat, vc_hat since they don't maintain rigid
        # body constraint yet:
        for k in range(ncorners):
            xk_n, yk_n, uk_n, vk_n = corners_n[k,:]  # unpack k'th row
            xc_hat.append(xk_n + dt*uk_n)
            yc_hat.append(yk_n + dt*vk_n)

        # remap to original shape, maintaining rigidity:
        xc_np1, yc_np1, theta_np1 = remap(xc_hat, yc_hat, z_np1, 
                                          debris.get_corners)
        # corresponding z vector for new position:
        z_np1 = (xc_np1[0], yc_np1[0], theta_np1)
    
        # compute new velocity at n+1 after remapping:
        
        # average fluid depth:

        hc_np1 = array([h(x,y,t_np1) for x,y in zip(xc_np1,yc_np1)])
        h_ave = hc_np1.mean()
        
        
        if h_ave < debris.draft:
            if (abs(corners_n[:,2]).max() + abs(corners_n[:,3]).max()) < 1e-3:
                # corner velocities at t_n were all zero, check static friction:
                friction = 'static'
            else:
                friction = 'kinetic'
            if verbose:
                print('%s friction at t = %.2f with h_ave = %.1f' \
                    % (friction,t_np1,h_ave))
        else:
            # debris is not touching ground:
            friction = None
            if verbose:
                print('No friction at t = %.2f with h_ave = %.1f' \
                    % (t_np1,h_ave))
        
        info_np1 = {'friction': friction}  # pass back for plotting purposes

        #print('At t = %.2f with h_ave = %.1f' % (t_np1,h_ave))
        wet_face_height = min(h_ave, debris.draft)
        face_area = debris.face_width * wet_face_height
        
        # split up area and mass between corners so each can be
        # accelerated separately
        corner_face_area = face_area / ncorners
        corner_bottom_area = debris.bottom_area / ncorners
        corner_mass = debris.mass / ncorners
        
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
                if friction is None:
                    # to advect with flow, corner velocity = fluid velocity:
                    uk_np1 = uk_f
                    vk_np1 = vk_f
                else:
                    # grounded:
                    uk_np1 = 0.
                    vk_np1 = 0.
            else:
                # compute forces and acceleration for this corner:
                
                # force Ffluid exerted by fluid based on velocity difference:
                du = uk_f - uk_n
                dv = vk_f - vk_n
                sk_f = sqrt(du**2 + dv**2)
                Ffluid_x = 0.5 * debris.rho_water * sk_f * du * corner_face_area
                Ffluid_y = 0.5 * debris.rho_water * sk_f * dv * corner_face_area
                #Ffluid_x = 0.
                #Ffluid_y = 0.
                Ffluid = sqrt(Ffluid_x**2 + Ffluid_y**2)
                
                if friction is not None:
                    Ffriction1 = debris.grav * corner_bottom_area \
                                * (debris.rho * debris.height - \
                                   debris.rho_water * h_ave)
                    if friction == 'static':
                        Ffriction = debris.friction_static * Ffriction1
                        Fnet = max(0., Ffluid - Ffriction)
                        if abs(Ffluid) < 0.01:
                            Fnet_x = Fnet_y = 0.
                        else:
                            Fnet_x = Ffluid_x * Fnet / Ffluid
                            Fnet_y = Ffluid_y * Fnet / Ffluid
                    elif friction == 'kinetic':
                        Ffriction = debris.friction_kinetic * Ffriction1
                        sk_n = sqrt(uk_n**2 + vk_n**2)
                        #Fnet_x = Ffluid_x - Ffriction * uk_n / sk_n
                        #Fnet_y = Ffluid_y - Ffriction * vk_n / sk_n
                        
                        Fnet_x = Ffluid_x
                        Fnet_y = Ffluid_y
                        
                    if verbose:
                        print('k = %i, Ffluid = %.3f, Ffriction = %.3f' \
                            % (k,Ffluid,Ffriction))

                else:
                    # not in contact with bottom, only fluid force:
                    Fnet_x = Ffluid_x
                    Fnet_y = Ffluid_y
                    
                uk_np1 = uk_n + dt*Fnet_x / corner_mass
                vk_np1 = vk_n + dt*Fnet_y / corner_mass
                
                if friction == 'kinetic':
                     decay = exp(- Ffriction / sk_n  * dt/corner_mass)
                     uk_np1 *= decay
                     vk_np1 *= decay
                
                if verbose:
                    print('k = %i, Fnet_x = %.3f, Fnet_y = %.3f' \
                        % (k,Fnet_x,Fnet_y))
            
            uc_np1[k] = uk_np1
            vc_np1[k] = vk_np1
        
        corners_np1 = vstack([xc_np1,yc_np1,uc_np1,vc_np1]).T
        corner_paths.append([t_np1, corners_np1, info_np1])
    
    return corner_paths         
    
def test_corner_paths_accel():
    
    debris = DebrisObject()
    debris.L = [1,1,1,1]
    debris.phi = [pi/2, pi/2, pi/2, pi/2]
    #debris.z0 = [3,1,pi/4]
    debris.z0 = [0,0,0]
    #debris.advect = False
    debris.rho = 900.
        
    u,v = velocities_shear()
    #h = lambda x,y,t: max(0., 0.2*(30.-t))  # fluid depth
    h = lambda x,y,t: where(t<30, 0.5*(1 - cos(2*pi*t/15.)), 0.)

    
    t0 = 0.
    nsteps = 50
    dt = 1.
    corner_paths_a = make_corner_paths_accel(debris,h,u,v,t0,dt,nsteps)
    corner_paths = make_corner_paths(debris,u,v,t0,dt,nsteps)
    
    figure(1);clf();
    for k in range(nsteps):
        tk,cpk,info = corner_paths[k]
        plot(cpk[:,0],cpk[:,1],'b')
        if mod(k,5) == 0:
            text(cpk[:,0].mean(), cpk[:,1].mean()+1, 't = %.1f' % tk,
                 color='b',fontsize=8)
        # one edge different color to show orientation:
        #plot(cpk[:2,0],cpk[:2,1],'c')
        tk,cpk,info = corner_paths_a[k]
        plot(cpk[:,0],cpk[:,1],'r')
        if mod(k,5) == 0:
            text(cpk[:,0].mean(), cpk[:,1].mean()-1.5, 't = %.1f' % tk,
                 color='r',fontsize=8)

    axis('equal')
    grid(True)
    
    return corner_paths, corner_paths_a


if __name__ == '__main__':

    debris = DebrisObject()
    debris.L = [1,1,1,1]
    debris.phi = [pi/2, pi/2, pi/2, pi/2]
    #debris.z0 = [3,1,pi/4]
    debris.z0 = [0,0,0]
    debris.advect = False
    debris.rho = 100.
    print('Draft = %.2fm' % debris.draft)
        
    u,v = velocities_shear()
    #h = lambda x,y,t: max(0., 0.2*(30.-t))  # fluid depth
    h = lambda x,y,t: where(t<300, 0.5*(1 - cos(2*pi*t/15.)), 0.)
    #h = lambda x,y,t: 0.55
    
    t0 = 0.
    nsteps = 150
    dt = 0.3
    corner_paths_a = make_corner_paths_accel(debris,h,u,v,t0,dt,nsteps)
    
    if 1:
        # Second object:
        debris2 = DebrisObject()
        debris2.L = 8 * [0.5]
        debris2.phi = [pi/2, 0., pi/2, 0., pi/2, 0., pi/2, 0.]
        debris2.z0 = [0,0,0]
        debris2.advect = False
        debris2.rho = 100.
        corner_paths_2 = make_corner_paths_accel(debris2,h,u,v,t0,dt,nsteps)
    else:
        corner_paths_2 = None
    
    fig = figure(1, figsize=(12,6))
    clf()
    tk,cpk,info = corner_paths_a[0]
    c = {None:'g', 'static':'r', 'kinetic':'orange'}
    plotk_a, = plot(cpk[:,0],cpk[:,1],color=c[info['friction']],lw=2)
    if corner_paths_2:
        tk,cpk2,info = corner_paths_2[0]
        plotk_2, = plot(cpk2[:,0],cpk2[:,1],color=c[info['friction']],lw=3)
    h0 = h(0,0,tk)
    title_text = title('time t = %.2fs, h = %.2fm' % (tk,h0))
    axis('scaled')
    axis([-1,20,-2,2])
    grid(True)
    
    def update(k):
        tk,cpk,info = corner_paths_a[k]
        plotk_a.set_data(cpk[:,0],cpk[:,1])
        plotk_a.set_color(c[info['friction']])
        if corner_paths_2:
            tk,cpk2,info = corner_paths_2[k]
            plotk_2.set_data(cpk2[:,0],cpk2[:,1])
            plotk_2.set_color(c[info['friction']])
        h0 = h(0,0,tk)
        title_text.set_text('time t = %.2fs, h = %.2fm' % (tk,h0))
        
    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=len(corner_paths_a), 
                                   interval=200, blit=False)

    fname_mp4 = 'corner_paths.mp4'
    fps = 5
    print('Making mp4...')
    animation_tools.make_mp4(anim, fname_mp4, fps)
    
