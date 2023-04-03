
import sys
if 'matplotlib' not in sys.modules:
    # for running on the cluster
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend


from pylab import *
import os
from matplotlib import colors

from clawpack.geoclaw.data import Rearth,DEG2RAD

    
coordinate_system = 1  # need to generalize

# =============================================
# Functions for tracking debris
# =============================================

# The functions make_fgout_fcn_xy and make_fgout_fcn_xyt are used below
# for interpolating the fluid velocity u,v from fgout frames to debris locations.
    

def make_fgout_fcn_xy(fgout, qoi, method='nearest',
                       bounds_error=False, fill_value=np.nan):
    """
    Create a function that can be called at (x,y) and return the qoi
    interpolated in space from the fgout array.
    
    qoi should be a string (e.g. 'u' or 'v') corresponding to 
    an attribute of fgout.
    
    The function returned takes arguments x,y that can be floats or 
    (equal length) 1D arrays of values that lie within the spatial
    extent of fgout. 
    
    bounds_error and fill_value determine the behavior if (x,y) is not in 
    the bounds of the data.
    """
    
    from scipy.interpolate import RegularGridInterpolator
    
    try:
        q = getattr(fgout,qoi)
    except:
        print('*** fgout missing attribute qoi = %s?' % qoi)
        
    err_msg = '*** q must have same shape as fgout.X\n' \
            + 'fgout.X.shape = %s,   q.shape = %s' % (fgout.X.shape,q.shape)
    assert fgout.X.shape == q.shape, err_msg
    
    x1 = fgout.X[:,0]
    y1 = fgout.Y[0,:]
    fgout_fcn1 = RegularGridInterpolator((x1,y1), q, method=method,
                bounds_error=bounds_error, fill_value=fill_value)

    def fgout_fcn(x,y):
        # evaluate at a single point or x,y arrays:
        xa = array(x)
        ya = array(y)
        xyout = vstack((xa,ya)).T
        qout = fgout_fcn1(xyout)
        if len(qout) == 1:
            qout = qout[0]  # return scalar
        return qout

    return fgout_fcn



def make_fgout_fcn_xyt(fgout1, fgout2, qoi, method_xy='nearest',
                       method_t='linear', bounds_error=False,
                       fill_value=np.nan):
    """
    Create a function that can be called at (x,y,t) and return the qoi
    interpolated in space and time between the two frames fgout1 and fgout2.
    
    qoi should be a string (e.g. 'u' or 'v') corresponding to 
    an attribute of fgout.
    
    method_xy is the method used in creating the spatial interpolator,
    method_t is the method used for interpolation in time.
    
    bounds_error and fill_value determine the behavior if (x,y,t) is not in 
    the bounds of the data.
    
    The function returned takes arguments x,y (floats or equal-length 1D arrays)
    of values that lie within the spatial extent of fgout1, fgout2
    (which are assumed to cover the same uniform grid at different times)
    and t should be a float that lies between fgout1.t and fgout2.t. 
    By default, the function returns np.nan outside of these limits.
    """
    
    assert allclose(fgout1.X, fgout2.X), '*** fgout1 and fgout2 must have same X'
    assert allclose(fgout1.Y, fgout2.Y), '*** fgout1 and fgout2 must have same Y'
    
    t1 = fgout1.t
    t2 = fgout2.t
    #assert t1 < t2, '*** expected fgout1.t < fgout2.t'
    
    fgout1_fcn_xy = make_fgout_fcn_xy(fgout1, qoi, method=method_xy,
                       bounds_error=bounds_error, fill_value=fill_value)
    fgout2_fcn_xy = make_fgout_fcn_xy(fgout2, qoi, method=method_xy,
                       bounds_error=bounds_error, fill_value=fill_value)
                
    def fgout_fcn(x,y,t):
        # function to evaluate at a single point or x,y arrays, at time t:
        xa = array(x)
        ya = array(y)
        tol = 1e-6  # to make sure it works ok when called with t=t1 or t=t2
        if t1-tol <= t <= t2+tol:
            alpha = (t-t1)/(t2-t1)
        elif bounds_error:
            errmsg = '*** argument t=%g should be between t1=%g and t2=%g' \
                     % (t,t1,t2)
            raise InputError(errmsg)
        else:
            qout = fill_value * ones(xa.shape)
            return qout

        qout1 = fgout1_fcn_xy(x,y)
        qout2 = fgout2_fcn_xy(x,y)

        if t1-1e-12 <= t <= t2+1e-12:
            alpha = (t-t1)/(t2-t1)
            qout = (1-alpha)*qout1 + alpha*qout2
        else:
            raise('*** t1=%g, t=%g, t2=%g' % (t1,t,t2))
            
        return qout
        
    return fgout_fcn
    


def move_debris(dbnos, debris_paths, fgout1, fgout2, drag_factor=None, 
                 grounding_depth=0.):
    """
    For each dbno in dbnos: debris_paths[dbno] is a 2D array and it is assumed
    that the last row has the form [t1, x1, y1, u1, v1] with the location and
    velocity of this debris particle at time t1, which should equal fgout1.t.
    
    Compute the location and velocity of the debris particle at time t2 and 
    append [t2, x2, y2, u2, v2] to the bottom of the array debris_paths[dbno].
    
    Currently implemented using the 2-step explicit Runge-Kutta method:
    1. Interpolate fgout1.u and fgout1.v to (x1,y1) and move the particle
    over time dt/2 with this velocity, to obtain (xm,ym).
    2. Interpolate (u,v) in space and time to this midpoint location, and then
    use this velocity to move the particle over time dt from (x1,y1) to (x2,y2).
    """
    
    t1 = fgout1.t
    t2 = fgout2.t
    dt = t2 - t1
    
    print('Moving debris over time dt = %g' % dt)
    print('       from t1 = %s to t2 = %.2f' % (t1,t2))

    h_fcn = make_fgout_fcn_xyt(fgout1, fgout2, 'h')
    u_fcn = make_fgout_fcn_xyt(fgout1, fgout2, 'u')
    v_fcn = make_fgout_fcn_xyt(fgout1, fgout2, 'v')

    for dbno in dbnos:
        debris_path = debris_paths[dbno]
        t1,xd1,yd1,ud1,vd1 = debris_path[-1,:]
        errmsg = '*** For dbno = %i, expected t1 = %.3f to equal fgout1.t = %.3f' \
                % (dbno, t1, fgout1.t)
        assert t1 == fgout1.t, errmsg
        
        try:
            gd = grounding_depth[dbno]  # if different for each particle
        except:
            gd = grounding_depth  # assume it's a scalar, same for all debris
            
        try:
            df = drag_factor[dbno]  # if different for each particle
        except:
            df = drag_factor  # assume it's a scalar, same for all debris
        
        if coordinate_system == 2:
            # x,y in degrees, u,v in m/s
            # convert u,v to degrees/second:
            ud1 = ud1 / (Rearth*DEG2RAD * cos(DEG2RAD*yd1))
            vd1 = vd1 / (Rearth*DEG2RAD)
            
        # Half time step with old velocities:
        xdm = xd1 + 0.5*dt*ud1
        ydm = yd1 + 0.5*dt*vd1
        
        tm = t1 + 0.5*dt  # t at midpoint in time
        
        # depth and fluid velocity at midpoint in time tm:
        hm = h_fcn(xdm,ydm,tm)
        um = u_fcn(xdm,ydm,tm)
        vm = v_fcn(xdm,ydm,tm)
    
        if hm < gd:
            # particle is grounded so velocities set to 0:
            udm = 0.
            vdm = 0.        
        elif df is None:
            # no drag factor and debris velocity = fluid velocity:
            udm = um
            vdm = vm
        else:
            # debris velocity (ud,vd) relaxes toward fluid velocity (u,v):
            # solve 
            #   d/dt ud = C * (u-ud)
            #   d/dt vd = C * (v-vd)
            # where C = df*sqrt((u-ud)**2 + (v-vd)**2).
            # Approximate this by setting C to value at t1, so simple
            # exponential decay with rate C:
            u1 = u_fcn(xd1,yd1,t1)
            v1 = v_fcn(xd1,yd1,t1)
            C = df*sqrt((u1-ud1)**2 + (v1-vd1)**2)
            udm = um - exp(-C*dt/2)*(um-ud1)
            vdm = vm - exp(-C*dt/2)*(vm-vd1)            
            
        if coordinate_system == 2:
            udm = udm / (Rearth*DEG2RAD * cos(DEG2RAD*ydm))
            vdm = vdm / (Rearth*DEG2RAD)
        
        # Take full time step with mid-point debris velocities:
        xd2 = xd1 + dt*udm
        yd2 = yd1 + dt*vdm
        x1b = 35.54
        x2b = 36.14
        y1b = 1.22
        y2b = y1b + 0.6
        dt2 = dt
        ncut = 0
        while (xd2>=x1b) and (xd2<=x2b) and (yd2>=y1b) and (yd2<=y2b) and (ncut<10):
            dt2 = dt2/2
            xd2 = xd1 + dt2*udm
            yd2 = yd1 + dt2*vdm
            ncut += 1
        if ncut >= 10:
            print('*** failed to reduce dt by enough')

        
        # Depth and fluid velocity at final time t2:
        h2 = h_fcn(xd2,yd2,t2)
        u2 = u_fcn(xd2,yd2,t2)
        v2 = v_fcn(xd2,yd2,t2)

        if h2 < gd:
            # particle is grounded so velocities set to 0:
            ud2 = 0.
            vd2 = 0.        
        elif df is None:
            # no drag factor and debris velocity = fluid velocity:
            ud2 = u2
            vd2 = v2
        else:
            # debris velocity (ud,vd) relaxes toward fluid velocity (u,v).
            # Take another half time step of decay from (udm,vdm),
            # now approximating C at the midpoint in time:
            C = df*sqrt((um-udm)**2 + (vm-vdm)**2)
            ud2 = u2 - exp(-C*dt/2)*(u2-udm)
            vd2 = v2 - exp(-C*dt/2)*(v2-vdm)    
        
        debris_paths[dbno] = vstack((debris_path, array([t2,xd2,yd2,ud2,vd2])))
        
    return debris_paths


            
def make_debris_paths(fgout_grid, fgframes, debris_paths, dbnos,
                      drag_factor=None, grounding_depth=0.):
    """
    dbnos is a list of debris particle labels (integers) to operate on.
    debris_paths a dictionary indexed by integer dbno.
    Each element 
        debris_path = debris_paths[dbno] 
    is a 2D numpy array with at least one row and three columns t,x,y.
    The last row of debris_path defines the starting time and location of 
    the particle.
    This routine loops over all fgout frames in fgframes, reads in the frame
    for fgout grid number fgno as fgout2, 
    and then calls move_debris to move each particle from time
    fgout1.t to fgout2.t, where fgout1 is the previous frame.  The new time
    and location are added as a new row in the debris_path array.
    """
    
    fgout1 = fgout_grid.read_frame(fgframes[0])
    
    for fgframe in fgframes[1:]:
        print('Trying to read fgno=%i, fgframe=%i' % (fgout_grid.fgno,fgframe))
        try:
            fgout2 = fgout_grid.read_frame(fgframe)
        except:
            print('Could not read file, exiting loop')
            break
        debris_paths = move_debris(dbnos, debris_paths, fgout1, fgout2,
                                   drag_factor=drag_factor, 
                                   grounding_depth=grounding_depth)
        fgout1 = fgout2
    return debris_paths

        
def plot_debris_1(t, debris_paths, dbnos, ax, 
                color='k', marker='s', markersize=4):
    """
    Plot the location of each debris particle from the list dbnos at time t,
    assuming debris_paths[dbno] has a row corresponding to this time.
    This assumes debris_paths contains the full debris paths as computed by
    make_debris_paths.
    """
    
    for dbno in dbnos:
        db = debris_paths[dbno]
        try:
            j = where(abs(db[:,0]-t) < 1e-6)[0].max()
        except:
            print('Did not find path for dbno=%i at t = %.3f' % (dbno,t))
            j = -1
        if j > -1:
            x = db[j,1]
            y = db[j,2]
            points, = ax.plot([x], [y], color=color, marker=marker, 
                             markersize=markersize)
        return points
            
def plot_debris(t, debris_paths, dbnos, ax, 
                color='k', marker='s', markersize=4):
    """
    Plot the location of each debris particle from the list dbnos at time t,
    assuming debris_paths[dbno] has a row corresponding to this time.
    This assumes debris_paths contains the full debris paths as computed by
    make_debris_paths.
    """
    
    xd = []
    yd = []
    for dbno in dbnos:
        db = debris_paths[dbno]
        try:
            j = where(abs(db[:,0]-t) < 1e-6)[0].max()
        except:
            print('Did not find path for dbno=%i at t = %.3f' % (dbno,t))
            j = -1
        if j > -1:
            xd.append(db[j,1])
            yd.append(db[j,2])
    points, = ax.plot(xd, yd, color=color, linestyle='', marker=marker, 
                     markersize=markersize)
    return points

# ===========
# debris pairs


def move_debris_pairs(dbnosA, dbnosB, debris_paths, fgout1, fgout2,
                      drag_factor=None, grounding_depth=0.):
    """
    The lists dbnosA and dbnosB should be the same length and the two particles
    dbnosA[k] and dbnosB[k] should be constrained to maintain constant distance 
    between them.
    """
    
    from clawpack.geoclaw.util import haversine
    
    if coordinate_system == 1:
        distfunc = lambda x1,y1,x2,y2: sqrt((x1-x2)**2 + (y1-y2)**2)

    elif coordinate_system == 2:
        distfunc = lambda x1,y1,x2,y2: haversine(x1,x2,y1,y2)
        

    # first move each particle by the unconstrained algorithm:
    dbnosAB = list(dbnosA) + list(dbnosB)
    move_debris(dbnosAB, debris_paths, fgout1, fgout2,
                drag_factor=drag_factor, grounding_depth=grounding_depth)
        
    # pdb; pdb.set_trace()
    
    # constrain motion so that adjacent particles remain const distance apart:
    for k in range(len(dbnosA)):
        dbnoA = dbnosA[k]
        dbnoB = dbnosB[k]
        
        # previous positions before move_debris performed above:
        tA_old,xdA_old,ydA_old,udA_old,vdA_old = debris_paths[dbnoA][-2,:]
        tB_old,xdB_old,ydB_old,udB_old,vdB_old = debris_paths[dbnoB][-2,:]
        dist_old = distfunc(xdB_old, ydB_old, xdA_old, ydA_old)

        # new positions computed by move_debris:
        tA_new,xdA_new,ydA_new,udA_new,vdA_new = debris_paths[dbnoA][-1,:]
        tB_new,xdB_new,ydB_new,udB_new,vdB_new = debris_paths[dbnoB][-1,:]
        dist_new = distfunc(xdB_new, ydB_new, xdA_new, ydA_new)
        
        # now adjust so that new dist agrees with dist_old
        # keep midpoint xdm,ydm fixed and adjust distance to each end:

        ratio = dist_old/dist_new

        xdm = 0.5*(xdA_new + xdB_new)
        ydm = 0.5*(ydA_new + ydB_new)
        dx_new = xdB_new - xdA_new
        dy_new = ydB_new - ydA_new
        xdA_new = xdm - ratio*0.5*dx_new
        ydA_new = ydm - ratio*0.5*dy_new
        xdB_new = xdm + ratio*0.5*dx_new
        ydB_new = ydm + ratio*0.5*dy_new
        dist_adjusted = distfunc(xdB_new, ydB_new, xdA_new, ydA_new)
        print('+++ dbnoA=%i: distances %.1f, %.1f, %.1f' \
                % (dbnoA, dist_old, dist_new, dist_adjusted))
        
        # Should reset ud, vd also!?
        
        debris_paths[dbnoA][-1,:] = tA_new,xdA_new,ydA_new,udA_new,vdA_new
        debris_paths[dbnoB][-1,:] = tB_new,xdB_new,ydB_new,udB_new,vdB_new
        
    return debris_paths


def make_debris_paths_pairs(fgout_grid, fgframes, debris_paths,
                            dbnosA, dbnosB, drag_factor=None, grounding_depth=0.):
    """
    dbnosA,dbnosB are equal-length lists of debris particle labels (integers)
    to operate on, specifying ed points of long debris particles.
    """
    
    fgout1 = fgout_grid.read_frame(fgframes[0])
    fgno = fgout_grid.fgno
    
    for fgframe in fgframes[1:]:
        print('Trying to read fgno=%i, fgframe=%i' % (fgno,fgframe))
        try:
            fgout2 = fgout_grid.read_frame(fgframe)
        except:
            print('Could not read file, exiting loop')
            break
        debris_paths = move_debris_pairs(dbnosA, dbnosB, debris_paths, 
                                         fgout1, fgout2,
                                         drag_factor=drag_factor, 
                                         grounding_depth=grounding_depth)
        fgout1 = fgout2
    return debris_paths


def plot_debris_pairs(t, debris_paths, dbnosA, dbnosB, ax,
                      color='k', linewidth=2):
    """
    Plot the location of each debris particle pair connected by a line,
    from the lists dbnosA, dbnosB at time t,
    assuming debris_paths[dbno] has a row corresponding to this time.
    This assumes debris_paths contains the full debris paths as computed by
    make_debris_paths.
    """
    
    for k in range(len(dbnosA)):
        dbnoA = dbnosA[k]
        dbnoB = dbnosB[k]
        dbA = debris_paths[dbnoA]
        dbB = debris_paths[dbnoB]
        try:
            j = where(abs(dbA[:,0]-t) < 1e-6)[0].max()
        except:
            print('Did not find path for dbno=%i at t = %.3f' % (dbno,t))
            j = -1
        if j > -1:
            xA = dbA[j,1]
            yA = dbA[j,2]
            xB = dbB[j,1]
            yB = dbB[j,2]
            ax.plot([xA,xB], [yA,yB], color=color, linewidth=linewidth)
            

    
def move_debris_substeps(dbnos, debris_paths, fgout1, fgout2, drag_factor=None, 
                 grounding_depth=0., mass=1e9, dradius=None, 
                 Kspring=None, tether=None, nsubsteps=1):
    """
    For each dbno in dbnos: debris_paths[dbno] is a 2D array and it is assumed
    that the last row has the form [t1, x1, y1, u1, v1] with the location and
    velocity of this debris particle at time t1, which should equal fgout1.t.
    
    Compute the location and velocity of the debris particle at time t2 and 
    append [t2, x2, y2, u2, v2] to the bottom of the array debris_paths[dbno].
    
    Currently implemented using the 2-step explicit Runge-Kutta method:
    1. Interpolate fgout1.u and fgout1.v to (x1,y1) and move the particle
    over time dt/2 with this velocity, to obtain (xm,ym).
    2. Interpolate (u,v) in space and time to this midpoint location, and then
    use this velocity to move the particle over time dt from (x1,y1) to (x2,y2).
    """

    from clawpack.geoclaw.util import haversine
    
    if coordinate_system == 1:
        distfunc = lambda x1,y1,x2,y2: sqrt((x1-x2)**2 + (y1-y2)**2)

    elif coordinate_system == 2:
        distfunc = lambda x1,y1,x2,y2: haversine(x1,x2,y1,y2)
        
        
    t1full = fgout1.t
    t2full = fgout2.t
    dt = t2full - t1full
    
    dt_substep = dt / nsubsteps
    
    print('Moving debris over time dt = %g with %i substeps' % (dt,nsubsteps))
    print('       from t1 = %s to t2 = %.2f' % (t1full,t2full))
    
    h_fcn = make_fgout_fcn_xyt(fgout1, fgout2, 'h')
    u_fcn = make_fgout_fcn_xyt(fgout1, fgout2, 'u')
    v_fcn = make_fgout_fcn_xyt(fgout1, fgout2, 'v')

    for ns in range(nsubsteps):
        ts1 = t1full + ns*dt_substep
        ts2 = ts1 + dt_substep
        #print('    substep %i with dt = %.3f starting at ts = %.3f' \
        #        % (ns,dt_substep,ts1))
        for dbno in dbnos:
            debris_path = debris_paths[dbno]
            t1,xd1,yd1,ud1,vd1 = debris_path[-1,:]
            errmsg = '*** For dbno = %i, expected t1 = %.3f to equal ts1 = %.3f, diff = %g' \
                    % (dbno, t1, ts1, t1-ts1)
            assert abs(t1-ts1)<1e-12, errmsg
            
            try:
                gd = grounding_depth[dbno]  # if different for each particle
            except:
                gd = grounding_depth  # assume it's a scalar, same for all debris
                
            try:
                df = drag_factor[dbno]  # if different for each particle
            except:
                df = drag_factor  # assume it's a scalar, same for all debris
                
            try:
                massdb = mass[dbno]
            except:
                massdb = mass  # assume it's a scalar, same for all debris
            
            if coordinate_system == 2:
                # x,y in degrees, u,v in m/s
                # convert u,v to degrees/second:
                ud1 = ud1 / (Rearth*DEG2RAD * cos(DEG2RAD*yd1))
                vd1 = vd1 / (Rearth*DEG2RAD)
            
        
            h1 = h_fcn(xd1,yd1,t1)
            u1 = u_fcn(xd1,yd1,t1)
            v1 = v_fcn(xd1,yd1,t1)
            
            # compute force on debris
            fxd = df*(u1 - ud1)  
            fyd = df*(v1 - vd1)
            
            if (dradius is not None):
                # compute inter-particle forces:
                for dbnok in dbnos:
                    if dbnok != dbno:
                        debris_pathk = debris_paths[dbnok]
                        tk,xdk,ydk,udk,vdk = debris_pathk[-1,:]
                        # need to fix for lat-lon:
                        dist = distfunc(xd1,yd1,xdk,ydk)
                        diamjk = dradius[dbno] + dradius[dbnok]
                        try:
                            Dt,Kt = tether(dbno,dbnok)
                        except:
                            Dt,Kt = nan,0.
                        if Kt > 0:
                            fxd = fxd + Kt*(xd1-xdk)*(Dt-dist)/dist
                            fyd = fyd + Kt*(yd1-ydk)*(Dt-dist)/dist
                        elif dist < diamjk:
                            fxd = fxd + Kspring*(xd1-xdk)*(diamjk-dist)/dist
                            fyd = fyd + Kspring*(yd1-ydk)*(diamjk-dist)/dist
            
            # acceleration:
            axd = fxd / massdb
            ayd = fyd / massdb
        
            
            ud2 = ud1 + dt_substep * axd
            vd2 = vd1 + dt_substep * ayd
            
            #print('+++ dbno = %4i, u1 = %.2f, v1 = %.2f' % (dbno,u1,v1))
            #print('+++ dbno = %4i, axd = %.2f, ayd = %.2f' % (dbno,axd,ayd))
            #print('+++ dbno = %4i, ud1 = %.2f, vd1 = %.2f' % (dbno,ud1,vd1))
            #print('+++ dbno = %4i, ud2 = %.2f, vd2 = %.2f' % (dbno,ud2,vd2))
                                                            
            # Take full time step with debris velocities:
            xd2 = xd1 + dt_substep*0.5*(ud1+ud2)
            yd2 = yd1 + dt_substep*0.5*(vd1+vd2)
            
            # Depth and fluid velocity at final time ts2:
            h2 = h_fcn(xd2,yd2,ts2)
            u2 = u_fcn(xd2,yd2,ts2)
            v2 = v_fcn(xd2,yd2,ts2)

            if h2 < gd:
                # particle is grounded so velocities set to 0:
                ud2 = 0.
                vd2 = 0.        
            elif df is None:
                # no drag factor and debris velocity = fluid velocity:
                ud2 = u2
                vd2 = v2
            
            #print('+++ saving for dbno=%i at ts2 = %.3f' % (dbno,ts2))
            debris_paths[dbno] = vstack((debris_path, array([ts2,xd2,yd2,ud2,vd2])))
        
    return debris_paths


            
def make_debris_paths_substeps(fgout_grid, fgframes, debris_paths, dbnos,
                      drag_factor=None, grounding_depth=0.,
                      mass=1e9, dradius=None, 
                      Kspring=None, tether=None, nsubsteps=1):
    """
    dbnos is a list of debris particle labels (integers) to operate on.
    debris_paths a dictionary indexed by integer dbno.
    Each element 
        debris_path = debris_paths[dbno] 
    is a 2D numpy array with at least one row and three columns t,x,y.
    The last row of debris_path defines the starting time and location of 
    the particle.
    This routine loops over all fgout frames in fgframes, reads in the frame
    for fgout grid number fgno as fgout2, 
    and then calls move_debris to move each particle from time
    fgout1.t to fgout2.t, where fgout1 is the previous frame.  The new time
    and location are added as a new row in the debris_path array.
    """
    
    fgout1 = fgout_grid.read_frame(fgframes[0])
    
    for fgframe in fgframes[1:]:
        print('Trying to read fgno=%i, fgframe=%i' % (fgout_grid.fgno,fgframe))
        try:
            fgout2 = fgout_grid.read_frame(fgframe)
        except:
            print('Could not read file, exiting loop')
            break
        debris_paths = move_debris_substeps(dbnos, debris_paths, fgout1, fgout2,
                                            drag_factor, grounding_depth, mass, 
                                            dradius, Kspring, tether, nsubsteps)
        fgout1 = fgout2
    return debris_paths
