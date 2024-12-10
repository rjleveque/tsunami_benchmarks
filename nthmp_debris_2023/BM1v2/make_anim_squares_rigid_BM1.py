
from pylab import *
import matplotlib.animation as animation
from matplotlib import colors
import numpy as np
from clawpack.visclaw import animation_tools, plottools, geoplot
from clawpack.geoclaw import fgout_tools
import debris_tracking



debris = debris_tracking.DebrisObject()
debris.L = 4 * [0.6]
debris.phi = [0, pi/2, pi/2, pi/2]
debris.height = 0.4
debris.bottom_area = debris.L[0]*debris.L[1]  # assuming rectangle
debris.face_width = debris.L[0]  # assuming square
debris.z0 = [34.34, 0.22, 0]
debris.friction_static = 0. #0.35
debris.friction_kinetic = 0. #0.25
debris.advect = True  #False
mass = 14.5 # kg
debris.rho = mass / (debris.height * debris.bottom_area)
print('Draft = %.2fm' % debris.draft)

use_sim_data = True

if not use_sim_data:

    from load_fgout import fgout_grid, h_fcn, u_fcn, v_fcn, fgout_times, \
                           fgout_h_txy, fgout_u_txy, fgout_v_txy
    fgout_grid_extent = fgout_grid.extent_edges    
               
else:
    # not working: u_vel_fcn((34.1, 34.34)) is nan
    # use provided sim data instead of GeoClaw fgout results:
    import pickle
    with open('../BM1/sim_data.pickle','rb') as f:
        sim_data = pickle.load(f)
    zeta_fcn = sim_data['zeta_fcn']  # function that interpolates zeta to (t,x)
    u_vel_fcn = sim_data['u_vel_fcn']  # function that interpolates u_vel to (t,x)

    def u_fcn(x,y,t):
        tx = (t,x)
        u_vel = float(u_vel_fcn(tx))
        if isnan(u_vel):
            u_vel = 0.
        return u_vel
        
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
        tx = (t,x)
        zeta = float(zeta_fcn(tx))
        if isnan(zeta):
            h = 0.
        else:
            h = max(zeta - B(x), 0.)
        return h
        
    fgout_times = array([0, 30., 60.])
    fgout_h_txy = zeros((3,400,60))
    fgout_u_txy = zeros((3,400,60))
    fgout_v_txy = zeros((3,400,60))
    fgout_grid_extent = [33.75, 43.75, -3.0, 3.0]
    



t0 = 34.
nsteps = 120 #221
dt = 0.1
corner_paths = debris_tracking.make_corner_paths_accel(debris,h_fcn,u_fcn,v_fcn,
                                                   t0,dt,nsteps,verbose=True)

corner_paths_2 = None  # no comparison

if 1:
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
        
    else:
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

    # plot debris square:
    
    c = {'no':'g', 'static':'r', 'kinetic':'orange'}
    
    tk,cpk,info = corner_paths[0]
    plotk, = plot(cpk[:,1],cpk[:,0],color=c[info['friction']],lw=2)
    if corner_paths_2:
        tk,cpk2,info = corner_paths_2[0]
        plotk_2, = plot(cpk2[:,1],cpk2[:,0],color=c[info['friction']],lw=3)

    
    def update(k):
            
        tk,cpk,info = corner_paths[k]
        plotk.set_data(cpk[:,1],cpk[:,0])
        plotk.set_color(c[info['friction']])
        if corner_paths_2:
            tk,cpk2,info = corner_paths_2[k]
            plotk_2.set_data(cpk2[:,1],cpk2[:,0])
            plotk_2.set_color(c[info['friction']])

        # color image:
        # choose closes fgout frame for now...
        fgframeno = where(fgout_times <= tk)[0].max()
        
        if imqoi == 'Depth':
            fgout_h = fgout_h_txy[fgframeno,:,:]
            eta_water = np.ma.masked_where(fgout_h < 1e-3, fgout_h)
            im.set_data(flipud(eta_water))
        else:
            fgout_s = sqrt(fgout_u_txy[fgframeno,:,:]**2 + \
                           fgout_v_txy[fgframeno,:,:]**2)
            im.set_data(flipud(fgout_s))
            
        t_str = '%.2f seconds' % tk
        title_text.set_text('%s at t = %s' % (imqoi,t_str))
        

    if 0:
        # plot center of mass path
        xcm = zeros(times.shape)
        ycm = zeros(times.shape)
        ndebris = 4
        for dbno in dbnosA:
            if dbno > 90:
                continue
            xdb = zeros(times.shape)
            ydb = zeros(times.shape)
            for corner in range(0,4):
                dbnoc = dbno + corner*1000
                xdb += debris_paths[dbnoc][:,1]
                ydb += debris_paths[dbnoc][:,2]
            xdb = xdb/4.
            ydb = ydb/4.
            xcm += minimum(xdb,43.75)
            ycm += ydb
            #import pdb; pdb.set_trace()
            ax.plot(ydb,xdb,'--',color='g',linewidth=0.5)
            
            dbxyt = vstack((times,xdb,ydb)).T
            fname = 'db%sxyt.txt' % str(dbno).zfill(2)
            savetxt(fname,dbxyt)
            print('Created ',fname)
            
        xcm = xcm/ndebris
        ycm = ycm/ndebris
        ax.plot(ycm,xcm,'-',color='r',linewidth=0.5)
        
        cmxyt = vstack((times,xcm,ycm)).T
        fname = 'cmxyt.txt'
        savetxt(fname,cmxyt)
        print('Created ',fname)

    
def plot_centroids(corner_paths):
    centroids_x = array([c[1][:,0].mean() for c in corner_paths])
    centroids_y = array([c[1][:,1].mean() for c in corner_paths])
    centroids_u = array([c[1][:,2].mean() for c in corner_paths])
    centroids_v = array([c[1][:,3].mean() for c in corner_paths])
    times = array([c[0] for c in corner_paths])
    
    d = loadtxt('/Users/rjl/git/tsunami_benchmarks/nthmp_debris_2023/BM1/Benchmark_1/comparison_data/paths_and_velocities/config1_vel.txt')
    figure(201,figsize=(8,7)); clf()
    subplot(211)
    plot(times,centroids_x,'b',label='GeoClaw')
    plot(d[:,0], d[:,1], 'c', label='provided data')
    legend(loc='upper left', framealpha=1)
    grid(True)
    title('x-position')
    ylabel('x (m)')
    subplot(212)
    plot(times,centroids_u,'b',label='GeoClaw')
    plot(d[:,0], d[:,3], 'c', label='provided data')
    legend(loc='upper left', framealpha=1)
    grid(True)
    title('cross-shore velocity')
    xlabel('time (sec)')
    ylabel('u velocity (m/s)')
    tight_layout()
    fname = 'centroids_xu.png'
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)
    
    figure(202, figsize=(5,7)); clf()
    plot(centroids_y, centroids_x, 'b-')
    plot(centroids_y[::3], centroids_x[::3], 'b.', markersize=3)
    axis([-1,1,44,34])
    title('Centroid location for config 1')
    xlabel('y (m)')
    ylabel('x (m)')
    grid(True)
    fname = 'centroids_yx.png'
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)

    return centroids_x, centroids_y, centroids_u, centroids_v
    
if __name__ == '__main__':
        
    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=len(corner_paths), 
                                   interval=200, blit=False)

    fname_mp4 = 'config1a.mp4'
    fps = 5
    print('Making mp4...')
    animation_tools.make_mp4(anim, fname_mp4, fps)
    
    centroids = plot_centroids(corner_paths)
    
    if use_sim_data:
        print('USING SIM DATA')
