
from pylab import *
import matplotlib.animation as animation
from matplotlib import colors
import numpy as np
from clawpack.visclaw import animation_tools, plottools, geoplot
from clawpack.geoclaw import fgout_tools
import debris_tracking
from shapely import Polygon
import copy, os

xplot = linspace(33.7625, 43.7375, 400)
yplot = linspace(-2.95, 2.95, 60)
Xplot,Yplot = meshgrid(xplot, yplot, indexing='ij')
fgout_grid_extent = [33.75,43.75,-3,3]

def track_debris(params):
    config = params['config']
    outdir = params['outdir']

    # small debris:
    debris = debris_tracking.DebrisObject()
    debris.L = 3 * [0.6]  # 3 sides define square, will have 4 corners
    debris.phi = [0, pi/2, pi/2]
    debris.height = 0.4
    debris.bottom_area = debris.L[0]*debris.L[1]  # assuming rectangle
    debris.face_width = debris.L[0]  # assuming square
    debris.friction_static = 0.35
    debris.friction_kinetic = 0.25
    debris.info = {'friction':'static'}
    debris.advect = False
    mass = 14.5 # kg
    debris.rho = mass / (debris.height * debris.bottom_area)
    print('Draft = %.2fm' % debris.draft)

    #config = 12
    use_sim_data = (outdir == 'sim_data')

    #plotdir = 'BM1_config12_dtopo4_block'
    #plotdir = 'BM1_config12_simdata'
    #os.system('mkdir -p %s' % plotdir)

    if config in [1,2]:
        debris_list = [debris]
        z0 = [34.34, 0.22, 0]  # determines initial debris location
        z0_list = [z0]


    if config in [3,4]:
        debris_list = [debris]
        z0 = [34.94, 0.82, 0]
        z0_list = [z0]

        debris = copy.deepcopy(debris)
        z0 = [34.94, 0.22, 0]
        debris_list.append(debris)
        z0_list.append(z0)

        debris = copy.deepcopy(debris)
        z0 = [34.34, 0.82, 0]
        debris_list.append(debris)
        z0_list.append(z0)

        debris = copy.deepcopy(debris)
        z0 = [34.34, 0.22, 0]
        debris_list.append(debris)
        z0_list.append(z0)

    if config == 12:
        # large debris
        debris = debris_tracking.DebrisObject()
        debris.L = 3 * [1.2]  # 3 sides define square, will have 4 corners
        debris.phi = [0, pi/2, pi/2]
        debris.height = 0.4
        debris.bottom_area = debris.L[0]*debris.L[1]  # assuming rectangle
        debris.face_width = debris.L[0]  # assuming square
        debris.friction_static = 0.35
        debris.friction_kinetic = 0.25
        debris.info = {'friction':'static'}
        debris.advect = False
        mass = 14.5 * 4 # kg
        debris.rho = mass / (debris.height * debris.bottom_area)
        print('Draft = %.2fm' % debris.draft)
        debris_list = [debris]

        z0 = [34.34, 0.22, 0]  # determines initial debris location
        z0_list = [z0]

    obst_list = []

    if 0:
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
        # domain defines back wall:
        domain = [0, 43.75, -3, 3]

    if config not in [1,3,7]:
        # obstacle -- stationary block:
        x1b = 35.54
        x2b = 36.14
        y1b = 1.22
        y2b = 1.82
        obst_x = array([x1b,x1b,x2b,x2b])
        obst_y = array([y1b,y2b,y2b,y1b])
        obst_p = vstack((obst_x,obst_y)).T
        obst_polygon = Polygon(obst_p)
        xcentroid = obst_x.mean()
        ycentroid = obst_y.mean()
        radius = 0.5*sqrt(2)*(x2b-x1b) # for square
        obst = {}
        obst['polygon'] = obst_polygon
        obst['xcentroid'] = xcentroid
        obst['ycentroid'] = ycentroid
        obst['radius'] = radius
        obst_list.append(obst)



    if not use_sim_data:

        #from load_fgout import fgout_grid, h_fcn, u_fcn, v_fcn, fgout_times, \
        #                       fgout_h_txy, fgout_u_txy, fgout_v_txy
        import load_fgout_fcn
        #if config == 3:
        #    outdir = '../BM1/SWE2d/_output_dtopo4_noblock'
        #elif config in [4,12]:
        #    #outdir = '../BM1/SWE2d/_output_dtopo4_noblock'
        #    outdir = '../BM1/SWE2d/_output_dtopo4_block'
        #else:
        #    print('Need to set outdir for config %i' % config)
        #fgframes = range(1,45)
        fgframes = None
        fgno = 1
        fgout_data = load_fgout_fcn.load_fgout(outdir,fgno,fgframes)
        # unpack:
        fgout_grid, h_fcn, u_fcn, v_fcn, fgout_times, \
               fgout_h_txy, fgout_u_txy, fgout_v_txy = fgout_data
        #fgout_grid_extent = fgout_grid.extent_edges
        #Xplot = fgout_grid.X
        #Yplot = fgout_grid.Y
        #fgout_grid_extent = fgout_grid.extent_edges

    else:
        # now working in way that we can plot h or u too
        # use provided sim data instead of GeoClaw fgout results:


        import pickle
        with open('../BM1/sim_data.pickle','rb') as f:
            sim_data = pickle.load(f)
        zeta_fcn = sim_data['zeta_fcn']  # function that interpolates zeta to (t,x)
        u_vel_fcn = sim_data['u_vel_fcn']  # function that interpolates u_vel to (t,x)

        def u_fcn(x,y,t):
            #tx = vstack((t,x)).T
            tx = [(t,xj) for xj in x.ravel()]
            u_vel = u_vel_fcn(tx)
            u_vel = reshape(u_vel, x.shape)
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
            #tx = vstack((t,x)).T
            tx = [(t,xj) for xj in x.ravel()]
            zeta = zeta_fcn(tx)
            zeta = reshape(zeta, x.shape)
            h = where(isnan(zeta), 0., zeta - B(x))
            return h

        if 0:
            # not used:
            fgout_times = array([0, 30., 60.])
            fgout_h_txy = zeros((3,400,60))
            fgout_u_txy = zeros((3,400,60))
            fgout_v_txy = zeros((3,400,60))
            #fgout_grid_extent = [33.75, 43.75, -3.0, 3.0]



    t0 = 34.
    nsteps = 221 #120 #221
    dt = 0.1
    times = arange(t0, t0+(nsteps+0.5)*dt, dt)

    debris_path_list = debris_tracking.make_debris_path_list(debris_list,
                                        z0_list,obst_list,domain,t0,dt,nsteps,
                                        h_fcn,u_fcn,v_fcn, verbose=False)


    # ===========
    # plotting



    bgimage = None
    #plot_extent = [34, 43.75, -3, 3]
    xlimits = [-3, 3]
    ylimits = [43.75, 34]
    #flipud = lambda A: fliplr(flipud(A))  # x is vertical in plots
    flipud = lambda A: A

    color = 'k'
    linewidth = 1

    fge1 = fgout_grid_extent  # for imshow plots
    fgout_extent = [fge1[2],fge1[3],fge1[1],fge1[0]]

    fig,ax = subplots(figsize=(6,7))

    ax.set_xlim(xlimits)
    ax.set_ylim(ylimits)

    #ax.plot([y1b,y1b,y2b,y2b,y1b], [x1b,x2b,x2b,x1b,x1b], 'g')

    imqoi = 'Depth'
    #imqoi = 'Speed'
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

        #fgout_h = fgout_h_txy[fgframeno,:,:]
        fgout_h = h_fcn(Xplot,Yplot,t0)
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

    else:
        # speed
        s_units = 'm/s'
        #fgout_s = sqrt(fgout_u_txy[fgframeno,:,:]**2 + \
        #               fgout_v_txy[fgframeno,:,:]**2)
        fgout_s = sqrt(u_fcn(Xplot,Yplot,t0)**2 + \
                       v_fcn(Xplot,Yplot,t0)**2)
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

    c = {'no':'g', 'static':'r', 'kinetic':'orange'}

    plot_debris_list = []
    for debris_path in debris_path_list:
        t_n = t0
        z_n = debris_path.z_path[0]
        info_n = debris_path.info_path[0]
        xc,yc = debris.get_corners(z_n, close_poly=True)
        plotn, = plot(yc, xc, color=c[info_n['friction']], lw=3)
        plot_debris_list.append(plotn)

    #plot_obst_list = []
    for obst in obst_list:
        obst_xy = array(obst['polygon'].exterior.coords)
        plot(obst_xy[:,1], obst_xy[:,0], 'k')

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
            plotn.set_color(c[info_n['friction']])


        # color image:
        # choose closes fgout frame for now...
        #fgframeno = where(fgout_times <= t_n)[0].max()

        if imqoi == 'Depth':
            #fgout_h = fgout_h_txy[fgframeno,:,:]
            fgout_h = h_fcn(Xplot,Yplot,t_n)
            eta_water = np.ma.masked_where(fgout_h < 1e-3, fgout_h)
            im.set_data(flipud(eta_water))
        else:
            #fgout_s = sqrt(fgout_u_txy[fgframeno,:,:]**2 + \
            #               fgout_v_txy[fgframeno,:,:]**2)
            fgout_s = sqrt(u_fcn(Xplot,Yplot,t_n)**2 + \
                           v_fcn(Xplot,Yplot,t_n)**2)
            im.set_data(flipud(fgout_s))

        t_str = '%.2f seconds' % t_n
        title_text.set_text('%s at t = %s (n=%i)' % (imqoi,t_str,n))

    return fig,update,debris_path_list



def plot_centroids(debris_path_list,params):
    config = params['config']
    outdir = params['outdir']
    plotdir = params['plotdir']
    use_sim_data = outdir == 'sim_data'

    os.system('mkdir -p %s' % plotdir)

    # need to fix for multiple debris
    for debris_path in debris_path_list[:1]:
        centroids_x = array([xc.mean() for xc in debris_path.x_path])
        centroids_y = array([yc.mean() for yc in debris_path.y_path])
        centroids_u = array([uc.mean() for uc in debris_path.u_path])
        centroids_v = array([vc.mean() for vc in debris_path.v_path])
        times = debris_path.times

        if use_sim_data:
            fluid_label = 'provided flowfield'
        else:
            fluid_label = 'GeoClaw SWE flowfield'

        d = loadtxt('/Users/rjl/git/tsunami_benchmarks/nthmp_debris_2023/BM1/Benchmark_1/comparison_data/paths_and_velocities/config%i_vel.txt' \
            % config)
        figure(201,figsize=(8,8)); clf()
        xlimits = (30,55)

        subplot(311)
        plot(times,centroids_x,'b',label='GeoClaw tracking with %s' % fluid_label)
        plot(d[:,0], d[:,1], 'c', label='provided comparison data')
        legend(loc='lower right', framealpha=1)
        grid(True)
        title('x-position for config %i' % config)
        ylabel('x (m)')
        xlim(xlimits)
        ylim(34,44)

        subplot(312)
        plot(times,centroids_u,'b',label='GeoClaw tracking with %s' % fluid_label)
        plot(d[:,0], d[:,3], 'c', label='provided comparison data')
        legend(loc='upper right', framealpha=1)
        grid(True)
        title('cross-shore velocity')
        ylabel('u velocity (m/s)')
        xlim(xlimits)
        ylim(-1,2)

        if 0:
            subplot(313)
            centroids_h = h_fcn(centroids_x, centroids_y, times)
            plot(times,100*centroids_h,'r',label=fluid_label)
            legend(loc='upper right', framealpha=1)
            grid(True)
            title('water depth')
            ylabel('water depth (cm)')
            xlim(xlimits)

        xlabel('time (sec)')

        tight_layout()
        fname = '%s/centroids_xu.png' % plotdir
        savefig(fname, bbox_inches='tight')
        print('Created ',fname)

    figure(202, figsize=(5,7)); clf()
    for debris_path in debris_path_list:
        centroids_x = array([xc.mean() for xc in debris_path.x_path])
        centroids_y = array([yc.mean() for yc in debris_path.y_path])
        centroids_u = array([uc.mean() for uc in debris_path.u_path])
        centroids_v = array([vc.mean() for vc in debris_path.v_path])
        times = debris_path.times

        plot(centroids_y, centroids_x, 'b-')
        plot(centroids_y[::3], centroids_x[::3], 'b.', markersize=3)
    axis([-2,2,44,34])
    title('Centroid location for config %i' % config)
    xlabel('y (m)')
    ylabel('x (m)')
    grid(True)
    fname = '%s/centroids_yx.png' % plotdir
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)

    return centroids_x, centroids_y, centroids_u, centroids_v

def plot_frames(fig,update,params,n_list=None):
    config = params['config']
    outdir = params['outdir']
    plotdir = params['plotdir']
    use_sim_data = outdir == 'sim_data'

    os.system('mkdir -p %s' % plotdir)
    if n_list is None:
        if config > 0:
            n_list = range(0,221,20)
    for n in n_list:
        update(n)
        fname = '%s/frame%s.png' % (plotdir,str(n).zfill(4))
        fig.savefig(fname, bbox_inches='tight')
        print('Saved ', fname)

def make_animation(fig,update,frames,params):
    config = params['config']
    plotdir = params['plotdir']

    os.system('mkdir -p %s' % plotdir)
    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=frames,
                                   interval=200, blit=False)

    fname_mp4 = '%s/BM1_config%s.mp4' % (plotdir,config)
    fps = 5
    print('Making mp4...')
    animation_tools.make_mp4(anim, fname_mp4, fps)


if __name__ == '__main__':

    if 0:
        print('Making anim...')
        anim = animation.FuncAnimation(fig, update,
                                       frames=len(debris_path.times),
                                       interval=200, blit=False)

        fname_mp4 = '%s/BM1_config%s.mp4' % (plotdir,config)
        fps = 5
        print('Making mp4...')
        animation_tools.make_mp4(anim, fname_mp4, fps)

        plot_frames()
        centroids = plot_centroids(debris_path_list)

        if use_sim_data:
            print('USING SIM DATA')

    run_all = False

    if 1 or run_all:
        params = {}
        params['config'] = 1
        params['outdir'] = 'sim_data'
        params['plotdir'] = 'BM1_config%i_simdata' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)

    if 0 or run_all:
        params = {}
        params['config'] = 1
        params['outdir'] = '../BM1/SWE2d/_output_noblock'
        params['plotdir'] = 'BM1_config%i_noblock' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)

    if 1 or run_all:
        params = {}
        params['config'] = 1
        params['outdir'] = '../BM1/SWE2d/_output_dtopo1'
        params['plotdir'] = 'BM1_config%i_dtopo1' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)


    if 0 or run_all:
        params = {}
        params['config'] = 3
        params['outdir'] = 'sim_data'
        params['plotdir'] = 'BM1_config%i_simdata' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)

    if 0 or run_all:
        params = {}
        params['config'] = 3
        params['outdir'] = '../BM1/SWE2d/_output_noblock'
        params['plotdir'] = 'BM1_config%i_noblock' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)

    if 0 or run_all:
        params = {}
        params['config'] = 3
        params['outdir'] = '../BM1/SWE2d/_output_dtopo4'
        params['plotdir'] = 'BM1_config%i_dtopo4' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)

    if 0 or run_all:
        params = {}
        params['config'] = 4
        params['outdir'] = 'sim_data'
        params['plotdir'] = 'BM1_config%i_simdata' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)

    if 0 or run_all:
        params = {}
        params['config'] = 4
        params['outdir'] = '../BM1/SWE2d/_output_block'
        params['plotdir'] = 'BM1_config%i_block' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)

    if 0 or run_all:
        params = {}
        params['config'] = 4
        params['outdir'] = '../BM1/SWE2d/_output_dtopo4_block'
        params['plotdir'] = 'BM1_config%i_dtopo4_block' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)


    if 0 or run_all:
        params = {}
        params['config'] = 12
        params['outdir'] = 'sim_data'
        params['plotdir'] = 'BM1_config%i_simdata' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)

    if 0 or run_all:
        params = {}
        params['config'] = 12
        params['outdir'] = '../BM1/SWE2d/_output_block'
        params['plotdir'] = 'BM1_config%i_block' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)

    if 0 or run_all:
        params = {}
        params['config'] = 12
        params['outdir'] = '../BM1/SWE2d/_output_dtopo4_block'
        params['plotdir'] = 'BM1_config%i_dtopo4_block' % params['config']
        fig,update,debris_path_list = track_debris(params)
        plot_frames(fig,update,params,n_list=None)
        make_animation(fig,update,None,params)
        plot_centroids(debris_path_list,params)
