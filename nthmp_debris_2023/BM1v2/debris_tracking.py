"""
Tools for debris tracking based on given velocity field.
Rigid body motion is supported.
Need to add interaction between debris Objects or with walls.
"""

from pylab import *


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

    def get_corners(self, z, close_poly=False):
        """
        convert z = [x0,y0,theta] into a list of corners, assuming (x0,y0) is
        starting corner and moving at angle theta (up from x-axis) to 2nd corner.
        Length of first side is L0, then turn through angle phi1, etc.

        If `close` is True, repeat the first point at the end for plotting
        purposes, but in general do not want to repeat a corner in computing
        forces or remapping to a rigid body so that distinct corners have equal
        weight.
        """
        x0,y0,theta = z # unpack
        ncorners = len(self.L)+1
        if close_poly:
            ncorners += 1
        xc = zeros(ncorners)
        yc = zeros(ncorners)
        xc[0] = x0
        yc[0] = y0
        phitot = theta
        #import pdb; pdb.set_trace()
        for k in range(len(self.L)):
            phitot = phitot + self.phi[k]
            xc[k+1] = xc[k] + self.L[k]*cos(phitot)
            yc[k+1] = yc[k] + self.L[k]*sin(phitot)
        if close_poly:
            xc[-1] = xc[0]
            yc[-1] = yc[0]

        return xc,yc


class DebrisPath():

    def __init__(self, debris=None, times=[0.], z0=(0,0,0)):
        self.debris = debris # for debris object being tracked
        self.times = times   # list of times
        self.z0 = z0
        if debris is not None:
            xcorners,ycorners = debris.get_corners(z0)
            ncorners = len(xcorners)
            self.x_path = nan*ones((len(times), ncorners))
            self.x_path[0,:] = xcorners
            self.y_path = nan*ones((len(times), ncorners))
            self.y_path[0,:] = ycorners
            self.u_path = zeros((len(times), ncorners))
            self.v_path = zeros((len(times), ncorners))
            self.z_path = nan*ones((len(times), 3))
            self.z_path[0,:] = z0
            self.info_path = []  # for list of dictionaries
        else:
            self.x_path = None  # array of shape (len(times), len(debris.corners))
            self.y_path = None
            self.u_path = None
            self.v_path = None
            self.z_path = None
            self.info_path = None

def make_debris_path(debris,z0,t0,dt,nsteps,h_fcn,u_fcn,v_fcn,verbose=False):

    # initial corner locations:
    xc,yc = debris.get_corners(z0)
    ncorners = len(xc)

    # initial velocities assumed to be zero:
    uc = zeros(ncorners)
    vc = zeros(ncorners)

    times = arange(t0, t0+(nsteps+0.5)*dt, dt)
    debris_path = DebrisPath(debris, times, z0)

    info0 = {'friction': 'static'}  # since not moving at t0
    debris_path.info_path = [info0]  # will append new dict at each time


    for n in range(nsteps):

        t_n = times[n]
        xc_n = debris_path.x_path[n,:]  # location at t_n
        yc_n = debris_path.y_path[n,:]
        uc_n = debris_path.u_path[n,:]  # velocity at t_n based on (delta x)/dt
        vc_n = debris_path.v_path[n,:]

        # compute new (provisional) velocities uc_hat, vc_hat based on forces
        # of fluid and friction. (May take xc,yc to points violating rigidity)

        # average fluid depth:

        hc_n = array([h_fcn(x,y,t_n) for x,y in zip(xc_n,yc_n)])
        h_ave = hc_n.mean()

        friction = 'no'
        if h_ave < debris.draft:
            if abs(uc_n).max() + abs(vc_n).max() < 1e-3:
                # corner velocities at t_n were all zero, check static friction:
                friction = 'static'
            else:
                friction = 'kinetic'

        if friction == 'static' and debris.friction_static == 0.:
            friction = 'no'
        if friction == 'kinetic' and debris.friction_kinetic == 0.:
            friction = 'no'
        if verbose:
            print('%s friction at t = %.2f with h_ave = %.1f' \
                % (friction,t_n,h_ave))

        info_np1 = {'friction': friction}  # pass back for plotting purposes

        #print('At t = %.2f with h_ave = %.1f' % (t_np1,h_ave))
        wet_face_height = min(h_ave, debris.draft)
        face_area = debris.face_width * wet_face_height

        # split up area and mass between corners so each can be
        # accelerated separately
        corner_face_area = face_area / ncorners
        corner_bottom_area = debris.bottom_area / ncorners
        corner_mass = debris.mass / ncorners

        #import pdb; pdb.set_trace()

        # initialize uc_hat,vc_hat vectors
        uc_hat = zeros(ncorners)
        vc_hat = zeros(ncorners)

        for k in range(ncorners):
            # loop over corners, setting uc_hat, vc_hat at each corner
            # based on fluid velocity and forces at that corner:

            xk_n = xc_n[k]
            yk_n = yc_n[k]
            uk_n = uc_n[k]
            vk_n = vc_n[k]

            # fluid velocities at this corner:
            uk_fluid = u_fcn(xk_n, yk_n, t_n)
            vk_fluid = v_fcn(xk_n, yk_n, t_n)

            #print('+++ k = %i, uk_f = %.3f  vk_f = %.3f' % (k,uk_f,vk_f))
            if isnan(uk_fluid):
                print('*** uk_fluid is nan at ', xk_n, yk_n, t_n)
                import pdb; pdb.set_trace()

            if debris.advect:
                if friction == 'no':
                    # to advect with flow, corner velocity = fluid velocity:
                    uk_hat = uk_fluid
                    vk_hat = vk_fluid
                else:
                    # grounded:
                    uk_hat = 0.
                    vk_hat = 0.
            else:
                # compute forces and acceleration for this corner:

                # force Ffluid exerted by fluid based on velocity difference:
                du = uk_fluid - uk_n
                dv = vk_fluid - vk_n
                sk_f = sqrt(du**2 + dv**2)
                Ffluid_x = 0.5 * debris.rho_water * sk_f * du * corner_face_area
                Ffluid_y = 0.5 * debris.rho_water * sk_f * dv * corner_face_area
                #Ffluid_x = 0.
                #Ffluid_y = 0.
                Ffluid = sqrt(Ffluid_x**2 + Ffluid_y**2)

                if friction != 'no':
                    Ffriction1 = debris.grav * corner_bottom_area \
                                * (debris.rho * debris.height - \
                                   debris.rho_water * h_ave)
                    if friction == 'static':
                        Ffriction = debris.friction_static * Ffriction1
                        Fnet = max(0., Ffluid - Ffriction)
                        # now rescale (Ffluid_x, Ffluid_y) vector
                        # to have length Fnet (net force after static friction)
                        if abs(Ffluid) < 0.01:
                            Fnet_x = Fnet_y = 0.
                        else:
                            Fnet_x = Ffluid_x * Fnet / Ffluid
                            Fnet_y = Ffluid_y * Fnet / Ffluid
                        print('+++s at t = %.1f, k = %i, Ffluid = %.3f, Ffriction = %.3f, Fnet_x = %.3f' \
                             % (t_n, k, Ffluid, Ffriction, Fnet_x))
                    elif friction == 'kinetic':
                        Ffriction = debris.friction_kinetic * Ffriction1
                        sk_n = sqrt(uk_n**2 + vk_n**2)
                        #Fnet_x = Ffluid_x - Ffriction * uk_n / sk_n
                        #Fnet_y = Ffluid_y - Ffriction * vk_n / sk_n

                        Fnet_x = Ffluid_x
                        Fnet_y = Ffluid_y
                        print('+++k at t = %.1f, k = %i, Ffluid = %.3f, Ffriction = %.3f, Fnet_x = %.3f' \
                             % (t_n, k, Ffluid, Ffriction, Fnet_x))

                    if verbose:
                        print('k = %i, Ffluid = %.3f, Ffriction = %.3f' \
                            % (k,Ffluid,Ffriction))

                else:
                    # not in contact with bottom, only fluid force:
                    Fnet_x = Ffluid_x
                    Fnet_y = Ffluid_y

                uk_hat = uk_n + dt*Fnet_x / corner_mass
                vk_hat = vk_n + dt*Fnet_y / corner_mass

                if friction == 'kinetic':
                     decay = exp(-Ffriction / sk_n  * dt/corner_mass)
                     uk_hat *= decay
                     vk_hat *= decay

                if verbose:
                    print('k = %i, Fnet_x = %.3f, Fnet_y = %.3f' \
                        % (k,Fnet_x,Fnet_y))

            uc_hat[k] = uk_hat
            vc_hat[k] = vk_hat

        # move corners based on this provisional velocity:
        xc_hat = xc_n + dt*uc_hat
        yc_hat = yc_n + dt*vc_hat

        # remap to original shape, maintaining rigidity:
        z_guess = debris_path.z_path[n]  # from previous step
        xc_np1, yc_np1, theta_np1 = remap(xc_hat, yc_hat, z_guess,
                                          debris.get_corners)

        # use actual distance moved to recompute velocities over last step:
        uc_np1 = (xc_np1 - xc_n) / dt
        vc_np1 = (yc_np1 - yc_n) / dt

        # check for bouncing off back wall -- NEED TO IMPROVE
        xwall2 = 43.75

        dx_wall = max(xc_np1) - xwall2
        if dx_wall > 0:
            xc_np1 = [x - dx_wall for x in xc_np1]
            uc_np1 = [-0.5*u for u in uc_np1]
            print('+++ negated uc at wall')
            #uc_np1 = [0. for u in uc_np1]

        # corresponding z vector for new position:
        z_np1 = (xc_np1[0], yc_np1[0], theta_np1)

        debris_path.x_path[n+1,:] = xc_np1
        debris_path.y_path[n+1,:] = yc_np1
        debris_path.u_path[n+1,:] = uc_np1
        debris_path.v_path[n+1,:] = vc_np1
        debris_path.z_path[n+1,:] = z_np1
        debris_path.info_path.append(info_np1)

    return debris_path

def make_debris_path_list(debris_list, obst_list, z0_list,
                          t0,dt,nsteps, h_fcn,u_fcn,v_fcn, verbose=False):

    debris_path_list = []
    for dbno in range(len(debris_list)):
        debris = debris_list[dbno]
        z0 = z0_list[dbno]
        # initial corner locations:
        xc,yc = debris.get_corners(z0)
        ncorners = len(xc)

        # initial velocities assumed to be zero:
        uc = zeros(ncorners)
        vc = zeros(ncorners)

        times = arange(t0, t0+(nsteps+0.5)*dt, dt)
        debris_path = DebrisPath(debris, times, z0)

        info0 = {'friction': 'static'}  # since not moving at t0
        debris_path.info_path = [info0]  # will append new dict at each time
        debris_path_list.append(debris_path)


    for n in range(nsteps):

        z_guess_list = []
        xc_hat_list = []
        yc_hat_list = []

        for dbno in range(len(debris_list)):
            debris = debris_list[dbno]
            debris_path = debris_path_list[dbno]

            # compute provisional  xc_hat, yc_hat
            t_n = times[n]
            xc_n = debris_path.x_path[n,:]  # location at t_n
            yc_n = debris_path.y_path[n,:]
            uc_n = debris_path.u_path[n,:]  # velocity at t_n based on (delta x)/dt
            vc_n = debris_path.v_path[n,:]

            # compute new (provisional) velocities uc_hat, vc_hat based on forces
            # of fluid and friction. (May take xc,yc to points violating rigidity)

            # average fluid depth:

            hc_n = array([h_fcn(x,y,t_n) for x,y in zip(xc_n,yc_n)])
            h_ave = hc_n.mean()

            friction = 'no'
            if h_ave < debris.draft:
                if abs(uc_n).max() + abs(vc_n).max() < 1e-3:
                    # corner velocities at t_n were all zero, check static friction:
                    friction = 'static'
                else:
                    friction = 'kinetic'

            if friction == 'static' and debris.friction_static == 0.:
                friction = 'no'
            if friction == 'kinetic' and debris.friction_kinetic == 0.:
                friction = 'no'
            if verbose:
                print('%s friction at t = %.2f with h_ave = %.1f' \
                    % (friction,t_n,h_ave))

            info_np1 = {'friction': friction}  # pass back for plotting purposes
            debris_path.info_path.append(info_np1)


            #print('At t = %.2f with h_ave = %.1f' % (t_np1,h_ave))
            wet_face_height = min(h_ave, debris.draft)
            face_area = debris.face_width * wet_face_height

            # split up area and mass between corners so each can be
            # accelerated separately
            corner_face_area = face_area / ncorners
            corner_bottom_area = debris.bottom_area / ncorners
            corner_mass = debris.mass / ncorners

            #import pdb; pdb.set_trace()

            # initialize uc_hat,vc_hat vectors
            uc_hat = zeros(ncorners)
            vc_hat = zeros(ncorners)

            for k in range(ncorners):
                # loop over corners, setting uc_hat, vc_hat at each corner
                # based on fluid velocity and forces at that corner:

                xk_n = xc_n[k]
                yk_n = yc_n[k]
                uk_n = uc_n[k]
                vk_n = vc_n[k]

                # fluid velocities at this corner:
                uk_fluid = u_fcn(xk_n, yk_n, t_n)
                vk_fluid = v_fcn(xk_n, yk_n, t_n)

                #print('+++ k = %i, uk_f = %.3f  vk_f = %.3f' % (k,uk_f,vk_f))
                if isnan(uk_fluid):
                    print('*** uk_fluid is nan at ', xk_n, yk_n, t_n)
                    import pdb; pdb.set_trace()

                if debris.advect:
                    if friction == 'no':
                        # to advect with flow, corner velocity = fluid velocity:
                        uk_hat = uk_fluid
                        vk_hat = vk_fluid
                    else:
                        # grounded:
                        uk_hat = 0.
                        vk_hat = 0.
                else:
                    # compute forces and acceleration for this corner:

                    # force Ffluid exerted by fluid based on velocity difference:
                    du = uk_fluid - uk_n
                    dv = vk_fluid - vk_n
                    sk_f = sqrt(du**2 + dv**2)
                    Ffluid_x = 0.5 * debris.rho_water * sk_f * du * corner_face_area
                    Ffluid_y = 0.5 * debris.rho_water * sk_f * dv * corner_face_area
                    #Ffluid_x = 0.
                    #Ffluid_y = 0.
                    Ffluid = sqrt(Ffluid_x**2 + Ffluid_y**2)

                    if friction != 'no':
                        Ffriction1 = debris.grav * corner_bottom_area \
                                    * (debris.rho * debris.height - \
                                       debris.rho_water * h_ave)
                        if friction == 'static':
                            Ffriction = debris.friction_static * Ffriction1
                            Fnet = max(0., Ffluid - Ffriction)
                            # now rescale (Ffluid_x, Ffluid_y) vector
                            # to have length Fnet (net force after static friction)
                            if abs(Ffluid) < 0.01:
                                Fnet_x = Fnet_y = 0.
                            else:
                                Fnet_x = Ffluid_x * Fnet / Ffluid
                                Fnet_y = Ffluid_y * Fnet / Ffluid
                            print('+++s at t = %.1f, k = %i, Ffluid = %.3f, Ffriction = %.3f, Fnet_x = %.3f' \
                                 % (t_n, k, Ffluid, Ffriction, Fnet_x))
                        elif friction == 'kinetic':
                            Ffriction = debris.friction_kinetic * Ffriction1
                            sk_n = sqrt(uk_n**2 + vk_n**2)
                            #Fnet_x = Ffluid_x - Ffriction * uk_n / sk_n
                            #Fnet_y = Ffluid_y - Ffriction * vk_n / sk_n

                            Fnet_x = Ffluid_x
                            Fnet_y = Ffluid_y
                            print('+++k at t = %.1f, k = %i, Ffluid = %.3f, Ffriction = %.3f, Fnet_x = %.3f' \
                                 % (t_n, k, Ffluid, Ffriction, Fnet_x))

                        if verbose:
                            print('k = %i, Ffluid = %.3f, Ffriction = %.3f' \
                                % (k,Ffluid,Ffriction))

                    else:
                        # not in contact with bottom, only fluid force:
                        Fnet_x = Ffluid_x
                        Fnet_y = Ffluid_y

                    uk_hat = uk_n + dt*Fnet_x / corner_mass
                    vk_hat = vk_n + dt*Fnet_y / corner_mass

                    if friction == 'kinetic':
                         decay = exp(-Ffriction / sk_n  * dt/corner_mass)
                         uk_hat *= decay
                         vk_hat *= decay

                    if verbose:
                        print('k = %i, Fnet_x = %.3f, Fnet_y = %.3f' \
                            % (k,Fnet_x,Fnet_y))

                uc_hat[k] = uk_hat
                vc_hat[k] = vk_hat

            # move corners based on this provisional velocity:
            xc_hat = xc_n + dt*uc_hat
            yc_hat = yc_n + dt*vc_hat

            xc_hat_list.append(xc_hat)
            yc_hat_list.append(yc_hat)
            z_guess_list.append(debris_path.z_path[n])


        # remap to original shape, maintaining rigidity and
        # avoiding collisions:
        xc_list, yc_list, z_list = remap_avoid(xc_hat_list, yc_hat_list,
                                          debris_list, obst_list, z_guess_list)

        for dbno in range(len(debris_list)):
            debris = debris_list[dbno]
            debris_path = debris_path_list[dbno]
            xc_n = debris_path.x_path[n,:]  # location at t_n
            yc_n = debris_path.y_path[n,:]
            xc_np1 = xc_list[dbno]  # new location from remap
            yc_np1 = yc_list[dbno]
            z_np1 = z_list[dbno]
            # use actual distance moved to recompute velocities over last step:
            uc_np1 = (xc_np1 - xc_n) / dt
            vc_np1 = (yc_np1 - yc_n) / dt

            # check for bouncing off back wall -- USE OBSTACLE

            debris_path.x_path[n+1,:] = xc_np1
            debris_path.y_path[n+1,:] = yc_np1
            debris_path.u_path[n+1,:] = uc_np1
            debris_path.v_path[n+1,:] = vc_np1
            debris_path.z_path[n+1,:] = z_np1

    return debris_path_list




def remap(xc_hat, yc_hat, z_guess, get_corners):
    """
    Adjust xc_hat, yc_hat to xc,yc so that original shape is preserved.
    z_guess is initial guess for z defining xc,yc.
    """

    from scipy.optimize import least_squares

    def F(z, *args, **kwargs):
        """
        Objective function for least squares fitting
        """
        from shapely import Polygon

        # compute current guess at corners from z:
        xc,yc = get_corners(z)

        # compute residuals between current corners and target:
        f = zeros(2*len(xc))
        for i in range(len(xc)):
            f[2*i] = xc[i] - xc_hat[i]
            f[2*i+1] = yc[i] - yc_hat[i]

        return f

    #print('+++ xc_hat: ',xc_hat)
    #print('+++ z_guess: ',z_guess)
    #print('+++ F: ',F(z_guess))

    result = least_squares(F,z_guess)
    xc,yc = get_corners(result['x'])
    theta = result['x'][2]
    #print('+++ remap result: ',result)
    return xc,yc,theta


def remap_avoid(xc_hat_list, yc_hat_list, debris_list, obst_list, z_guess_list):
    """
    Adjust xc_hat, yc_hat to xc,yc so that original shape is preserved.
    z_guess is initial guess for z defining xc,yc.

    Avoid collisions.
    """

    from scipy.optimize import least_squares

    z_guess_all = z_guess_list[0]
    for dbno in range(1, len(z_guess_list)):
        z_guess_all = z_guess_all + z_guess_list[dbno]

    def F(z_all, *args, **kwargs):
        """
        Objective function for least squares fitting
        """
        from shapely import Polygon

        f_total = array([], dtype=float)

        for dbno in range(len(debris_list)):
            debris = debris_list[dbno]
            z = z_all[3*dbno:3*dbno+3]
            xc_hat = xc_hat_list[dbno]
            yc_hat = yc_hat_list[dbno]

            # compute current guess at corners from z:
            xc,yc = debris.get_corners(z)

            # compute residuals between current corners and target:
            f = zeros(2*len(xc))
            for i in range(len(xc)):
                f[2*i] = xc[i] - xc_hat[i]
                f[2*i+1] = yc[i] - yc_hat[i]

            f_total = hstack((f_total, f))

            # debris-obstacle collisions:

            debris_polygon = Polygon(vstack((xc,yc)).T)

            for obno in range(len(obst_list)):
                obst_polygon = obst_list[obno]
                overlap_area = obst_polygon.intersection(debris_polygon).area
                f_total = hstack((f_total, 100*overlap_area))

            # debris-debris collisions:

            for dbno2 in range(dbno+1,len(debris_list)):
                debris2 = debris_list[dbno2]
                z2 = z_guess_list[dbno2]
                xc2,yc2 = debris.get_corners(z2)
                debris2_polygon = Polygon(vstack((xc2,yc2)).T)
                overlap_area = debris2_polygon.intersection(debris_polygon).area
                f_total = hstack((f_total, 100*overlap_area))
                #print('+++ dbno = %i, dbno2 = %i' % (dbno,dbno2))

        return f_total

    z_guess_all = z_guess_list[0]
    for dbno in range(1,len(z_guess_list)):
        z_guess_all = hstack((z_guess_all, z_guess_list[dbno]))

    #print('+++ xc_hat_list: ',xc_hat_list)
    #print('+++ z_guess_all: ',z_guess_all)
    #print('+++ F: ',F(z_guess_all))

    result = least_squares(F,z_guess_all)
    #print('+++ remap_avoid result: ',result)

    z_all = result['x']
    xc_list = []
    yc_list = []
    theta_list = []
    z_list = []
    for dbno in range(len(debris_list)):
        z = z_all[3*dbno:3*dbno+3]
        z_list.append(z)
        xc,yc = debris_list[dbno].get_corners(z)
        xc_list.append(xc)
        yc_list.append(yc)
        theta_list.append(z_all[3*dbno+2])
    return xc_list,yc_list,z_list


# ==========================================
# Testing with simple examples....  Need to clean up

def velocities_shear():
    u = lambda x,y,t: y
    v = lambda x,y,t: zeros(x.shape)
    return u,v

def test_debris_path():

    debris = DebrisObject()
    debris.L = [1,1,1]
    debris.phi = [pi/2, pi/2, pi/2]
    #debris.z0 = [3,1,pi/4]

    debris.advect = True
    debris.rho = 900.

    u,v = velocities_shear()
    #h = lambda x,y,t: max(0., 0.2*(30.-t))  # fluid depth
    #h = lambda x,y,t: where(t<30, 0.5*(1 - cos(2*pi*t/15.)), 0.)
    h = lambda x,y,t: 10.

    t0 = 0.
    nsteps = 10
    dt = 0.5
    #z0 = [0,-0.5,0]
    z0 = [0,0.5,0]
    debris_path = make_debris_path(debris,z0,t0,dt,nsteps,h,u,v)

    figure(1);clf();
    for n in range(nsteps+1):
        t_n = debris_path.times[n]
        z_n = debris_path.z_path[n]
        xc,yc = debris.get_corners(z_n, close_poly=True)
        plot(xc, yc, label='t = %.1f' % t_n)
        #if mod(k,5) == 0:
        if 0:
            text(debris_path.x_path[k].mean(),
                 debris_path.y_path[k].mean()+1, 't = %.1f' % t_n,
                 color='b',fontsize=8)

    #legend()
    axis('equal')
    grid(True)

    return debris_path

def test_debris_path_list():

    debris = DebrisObject()
    debris.L = [1,1,1]
    debris.phi = [pi/2, pi/2, pi/2]
    #debris.z0 = [3,1,pi/4]

    debris.advect = True
    debris.rho = 900.

    u,v = velocities_shear()
    #h = lambda x,y,t: max(0., 0.2*(30.-t))  # fluid depth
    #h = lambda x,y,t: where(t<30, 0.5*(1 - cos(2*pi*t/15.)), 0.)
    h = lambda x,y,t: 10.

    t0 = 0.
    nsteps = 10
    dt = 0.5
    #z0 = [0,-0.5,0]
    z0 = [0,0.5,0]

    debris_list = [debris]
    obst_list = []
    z0_list = [z0]

    debris_path_list = make_debris_path_list(debris_list, obst_list, z0_list,
                                            t0,dt,nsteps,h,u,v)

    debris_path = debris_path_list[0]

    figure(2);clf();
    for n in range(nsteps+1):
        t_n = debris_path.times[n]
        z_n = debris_path.z_path[n]
        xc,yc = debris.get_corners(z_n, close_poly=True)
        plot(xc, yc, label='t = %.1f' % t_n)
        #if mod(k,5) == 0:
        if 0:
            text(debris_path.x_path[k].mean(),
                 debris_path.y_path[k].mean()+1, 't = %.1f' % t_n,
                 color='b',fontsize=8)

    #legend()
    axis('equal')
    grid(True)

    return debris_path


if __name__ == '__main__':

    import matplotlib.animation as animation
    from clawpack.visclaw import animation_tools

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
