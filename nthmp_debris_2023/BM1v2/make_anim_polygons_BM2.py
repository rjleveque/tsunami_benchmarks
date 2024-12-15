from pylab import *
import matplotlib.animation as animation
from matplotlib import colors
import numpy as np
from numpy import random
from clawpack.visclaw import animation_tools, plottools, geoplot
from clawpack.geoclaw import fgout_tools
import debris_tracking
from shapely import Polygon
import copy

debris_wood = debris_tracking.DebrisObject()
debris_wood.L = 3 * [0.102]  # 3 sides define square, will have 4 corners
debris_wood.phi = [0, pi/2, pi/2]
debris_wood.height = 0.051
debris_wood.bottom_area = debris_wood.L[0]*debris_wood.L[1]  # assuming rectangle
debris_wood.face_width = debris_wood.L[0]  # assuming square
debris_wood.friction_static = 0.71  # wood
debris_wood.friction_kinetic = 0.5 * debris_wood.friction_static
debris_wood.advect = False
debris_wood.rho = 648. # wood
print('Draft of wood block = %.3fm' % debris_wood.draft)
debris_wood.info = {'material':'wood', 'friction':'static'}

debris_hdpe = copy.deepcopy(debris_wood)
debris_hdpe.friction_static = 0.38  # HDPE
debris_hdpe.friction_kinetic = 0.5 * debris_hdpe.friction_static
debris_hdpe.rho = 987. # HDPE
print('Draft of hdpe block = %.3fm' % debris_hdpe.draft)
debris_hdpe.info = {'material':'hdpe', 'friction':'static'}

debris_list = []
z0_list = []

x1 = 31.29
y1 = -2*0.153 - 0.051
randomize = True

nrows = 4
n_per_row = 5

if randomize:
    rr = random.uniform(size=nrows*n_per_row*3)
    fname = 'random.txt'
    savetxt(fname, rr)
    print('Created ',fname)

blockno = 0

for j in range(nrows):
    for i in range(n_per_row):
        if mod(blockno,2) == 0:
            debris = copy.deepcopy(debris_wood)
        else:
            debris = copy.deepcopy(debris_hdpe)

        # initial debris location and angle:
        x_ij = x1 + j*0.153
        y_ij = y1 + i*0.153
        theta_ij = 0.
        if randomize:
            x_ij = x_ij + .05*rr[3*blockno + 0]
            y_ij = y_ij + .05*rr[3*blockno + 1]
            theta_ij = rr[3*blockno + 2] * pi/2
        z0 = [x_ij,y_ij,theta_ij]  # determines initial debris location
        z0_list.append(z0)
        debris.z0 = z0
        debris_list.append(debris)
        blockno += 1

if randomize:
    obst_list = []
    xc_hat_list = []
    yc_hat_list = []
    z_guess_list = z0_list
    for debris in debris_list:
        xc_hat, yc_hat = debris.get_corners(debris.z0)
        xc_hat_list.append(xc_hat)
        yc_hat_list.append(yc_hat)
    xc_list,yc_list,z_list = debris_tracking.remap_avoid(xc_hat_list,
                             yc_hat_list, debris_list, obst_list, z_guess_list)
    z0_list = z_list

obst_list = []

if 1:
    # back wall:
    xwall2 = 43.75
    x1b = xwall2
    x2b = xwall2 + 1
    y1b = -5.
    y2b = 5.
    obst_x = array([x1b,x1b,x2b,x2b])
    obst_y = array([y1b,y2b,y2b,y1b])
    obst_p = vstack((obst_x,obst_y)).T
    obst_polygon = Polygon(obst_p)
    obst_list.append(obst_polygon)

if 1:
    # block 1:
    x1b = 35.29
    x2b = x1b + 0.4
    y1b = -0.6
    y2b = y1b + 0.4
    obst_x = array([x1b,x1b,x2b,x2b])
    obst_y = array([y1b,y2b,y2b,y1b])
    obst_p = vstack((obst_x,obst_y)).T
    obst_polygon = Polygon(obst_p)
    obst_list.append(obst_polygon)

    # block 2:
    x1b = 35.29
    x2b = x1b + 0.4
    y1b = 0.2
    y2b = y1b + 0.4
    obst_x = array([x1b,x1b,x2b,x2b])
    obst_y = array([y1b,y2b,y2b,y1b])
    obst_p = vstack((obst_x,obst_y)).T
    obst_polygon = Polygon(obst_p)
    obst_list.append(obst_polygon)

use_sim_data = True

if not use_sim_data:

    # NEED TO FIX FOR BM2
    from load_fgout import fgout_grid, h_fcn, u_fcn, v_fcn, fgout_times, \
                           fgout_h_txy, fgout_u_txy, fgout_v_txy
    fgout_grid_extent = fgout_grid.extent_edges

else:
    # not working: u_vel_fcn((34.1, 34.34)) is nan
    # use provided sim data instead of GeoClaw fgout results:
    import pickle
    with open('../BM2/sim_data.pickle','rb') as f:
        sim_data = pickle.load(f)
    zeta_fcn = sim_data['zeta_fcn']  # function that interpolates zeta to (t,x)
    u_vel_fcn = sim_data['u_vel_fcn']  # function that interpolates u_vel to (t,x)

    def u_fcn(x,y,t):
        tx = vstack((t,x)).T
        u_vel = u_vel_fcn(tx)
        u = where(isnan(u_vel), 0., u_vel)
        return u

    v_fcn = lambda x,y,t: 0.

    # piecewise linear bottom:
    #xB = array([0,10,17.5,32,43.75])
    #yB = array([0,0,0.5,1,1]) - 0.9017
    def B(x):
        x = array(x)
        B = zeros(x.shape)
        B = where(x<10, B, (x-10)/15.)
        B = where(x<17.5, B, 0.5 + (x-17.5)/29.)
        B = where(x<32, B, 1)
        B = B - 0.9017
        return B

    def h_fcn(x,y,t):
        tx = vstack((t,x)).T
        zeta = zeta_fcn(tx)
        h = where(isnan(zeta), 0., zeta - B(x))
        return h

    fgout_times = array([0, 30., 60.])
    fgout_h_txy = zeros((3,400,60))
    fgout_u_txy = zeros((3,400,60))
    fgout_v_txy = zeros((3,400,60))
    fgout_grid_extent = [33.75, 43.75, -3.0, 3.0]



t0 = 34.
nsteps = 71 #221
dt = 0.3

debris_path_list = debris_tracking.make_debris_path_list(debris_list,
                                    obst_list,z0_list,t0,dt,nsteps,
                                    h_fcn,u_fcn,v_fcn, verbose=False)



if 1:
    # ===========
    # plotting

    #fgout_grid_extent = [30, 43.75, -3, 3]
    fgout_grid_extent = [30, 40, -2, 2]

    bgimage = None
    #plot_extent = [34, 43.75, -3, 3]
    #flipud = lambda A: fliplr(flipud(A))  # x is vertical in plots
    flipud = lambda A: A

    color = 'k'
    linewidth = 1

    fge1 = fgout_grid_extent  # for imshow plots
    fgout_extent = [fge1[2],fge1[3],fge1[1],fge1[0]]
    xlimits = [fge1[2],fge1[3]]
    ylimits = [fge1[1],fge1[0]]

    fig,ax = subplots(figsize=(6,7))

    ax.set_xlim(xlimits)
    ax.set_ylim(ylimits)

    #ax.plot([y1b,y1b,y2b,y2b,y1b], [x1b,x2b,x2b,x1b,x1b], 'g')

    imqoi = 'Depth'
    fgframeno = 0  # initial fgout frame to use from fgout_h_txy

    if imqoi=='Depth':
        # depth

        a = 1.
        #cmap_depth = mpl.colors.ListedColormap([
        #                [.6,.6,1,a],[.3,.3,1,a],[0,0,1,a], [1,.8,.8,a],[1,.6,.6,a],
        #                [1,0,0,a]])
        cmap_depth = mpl.colors.ListedColormap([
                        [.8,.8,1,a],[.7,.7,1,a],[.6,.6,1,a],[.5,.5,1,a]])

        # Set color for value exceeding top of range:
        cmap_depth.set_over(color=[0.4,0.4,1,a])


        if bgimage:
            # Set color to transparent where s < 1e-3:
            cmap_depth.set_under(color=[1,1,1,0])
        else:
            # Set color to white where s < 1e-3:
            cmap_depth.set_under(color=[1,1,1,a])

        #bounds_depth = np.array([0,1,2,3,4,5])
        #bounds_depth = np.array([0,0.04,0.08,0.12,0.16,0.20,0.24])
        bounds_depth = np.array([0.001,0.04,0.06,0.08,0.10])

        norm_depth = colors.BoundaryNorm(bounds_depth, cmap_depth.N)

        fgout_h = fgout_h_txy[fgframeno,:,:]
        eta_water = np.ma.masked_where(fgout_h < 1e-3, fgout_h)

        im = imshow(flipud(eta_water), extent=fgout_extent,
                        #cmap=geoplot.tsunami_colormap)
                        cmap=cmap_depth, norm=norm_depth)
        im.set_clim(-5,5)

        cb = colorbar(im, extend='max', shrink=0.7)
        cb.set_label('meters')
        #contour(fgout.X, fgout.Y, fgout.B, [0], colors='g', linewidths=0.5)

        #ax.set_aspect(1./cos(ylat*pi/180.))
        ticklabel_format(useOffset=False)
        xticks(rotation=20)
        ax.set_xlim(xlimits)
        ax.set_ylim(ylimits)

        t_str = '%.2f seconds' % t0
        title_text = title('%s at t = %s' % (imqoi,t_str))

    elif imqoi=='Speed':
        # speed
        s_units = 'm/s'
        fgout_s = sqrt(fgout_u_txy[fgframeno,:,:]**2 + \
                       fgout_v_txy[fgframeno,:,:]**2)
        if s_units == 'knots':
            s = fgout_s * 1.9438  # convert m/s to knots
            bounds_speed = np.array([1e-3,3,4,6,9,12])  # knots
        else:
            s = fgout_s
            bounds_speed = np.array([1e-3,1.5,2,3,4.5,6])  # m/s
            bounds_speed = np.array([1e-3,0.1,0.2,0.3,0.4,0.5])  # m/s
            bounds_speed = np.array([1e-3,0.2,0.4,0.6,1.,1.5])  # m/s

        a = 1.
        cmap_speed = mpl.colors.ListedColormap([
                        [.3,.3,1,a],[0,0,1,a], [1,.8,.8,a],[1,.6,.6,a],
                        [1,0,0,a]])

        # Set color for value exceeding top of range to purple:
        cmap_speed.set_over(color=[1,0,1,a])

        if bgimage:
            # Set color to transparent where s < 1e-3:
            cmap_speed.set_under(color=[1,1,1,0])
        else:
            # Set color to white where s < 1e-3:
            cmap_speed.set_under(color=[1,1,1,a])

        norm_speed = colors.BoundaryNorm(bounds_speed, cmap_speed.N)

        im = imshow(flipud(s), extent=fgout_extent,
                    cmap=cmap_speed, norm=norm_speed)
        cb = colorbar(im, extend='max', shrink=0.7)
        cb.set_label(s_units)
        #contour(fgout.X, fgout.Y, fgout.B, [0], colors='g', linewidths=0.5)

        #ax.set_aspect(1./cos(ylat*pi/180.))
        ticklabel_format(useOffset=False)
        xticks(rotation=20)
        ax.set_xlim(xlimits)
        ax.set_ylim(ylimits)

        t_str = '%.2f seconds' % t0
        title_text = title('%s at t = %s' % (imqoi,t_str))

    # plot debris:

    #c = {'no':'g', 'static':'r', 'kinetic':'orange'}
    c = {'wood':'g', 'hdpe':'r'}

    plot_debris_list = []
    for debris_path in debris_path_list:
        t_n = t0
        z_n = debris_path.z_path[0]
        info_n = debris_path.info_path[0]
        xc,yc = debris.get_corners(z_n, close_poly=True)
        plotn, = plot(yc, xc, color=c[info_n['material']], lw=1)
        plot_debris_list.append(plotn)

    #plot_obst_list = []
    for obst_poly in obst_list:
        obst_xy = array(obst_poly.exterior.coords)
        plot(obst_xy[:,1], obst_xy[:,0], 'c')

    grid(True)

    def update(n):

        for dbno in range(len(debris_list)):
            debris_path = debris_path_list[dbno]
            plotn = plot_debris_list[dbno]
            t_n = debris_path.times[n]
            z_n = debris_path.z_path[n]
            info_n = debris_path.info_path[n]

            # update plot of debris:
            xc,yc = debris.get_corners(z_n, close_poly=True)
            plotn.set_data(yc,xc)
            plotn.set_color(c[info_n['material']])


        # color image:
        # choose closes fgout frame for now...
        fgframeno = where(fgout_times <= t_n)[0].max()

        if imqoi == 'Depth':
            fgout_h = fgout_h_txy[fgframeno,:,:]
            eta_water = np.ma.masked_where(fgout_h < 1e-3, fgout_h)
            im.set_data(flipud(eta_water))
        elif imqoi == 'Speed':
            fgout_s = sqrt(fgout_u_txy[fgframeno,:,:]**2 + \
                           fgout_v_txy[fgframeno,:,:]**2)
            im.set_data(flipud(fgout_s))

        t_str = '%.2f seconds' % t_n
        title_text.set_text('%s at t = %s' % (imqoi,t_str))


if __name__ == '__main__':

    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=len(debris_path.times),
                                   interval=200, blit=False)

    fname_mp4 = 'polygon_test_bm4.mp4'
    fps = 5
    print('Making mp4...')
    animation_tools.make_mp4(anim, fname_mp4, fps)

    #centroids = plot_centroids(debris_path_list)

    if use_sim_data:
        print('USING SIM DATA')
