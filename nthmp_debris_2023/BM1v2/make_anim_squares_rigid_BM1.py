
from pylab import *
import matplotlib.animation as animation
from matplotlib import colors
import numpy as np
from clawpack.visclaw import animation_tools, plottools, geoplot
from clawpack.geoclaw import fgout_tools
import RigidMotion
from load_fgout import h_fcn, u_fcn, v_fcn, fgout_times, fgout_h_txy

outdir = '../BM1/SWE2d/_output_dtopo1'
format = 'binary32'  # format of fgout grid output


print('Looking for output in ',outdir)

output_format = 'binary32'

# List of frames to use for making debris paths and animation:
#fgframes = range(10,121)
fgframes = range(10,21)
#fgframes = [30,31,32]


bgimage = None
#plot_extent = [34, 43.75, -3, 3]
xlimits = [-3, 3]
ylimits = [43.75, 34]
#flipud = lambda A: fliplr(flipud(A))  # x is vertical in plots
flipud = lambda A: A

color = 'k'
linewidth = 1


# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(1, outdir, output_format)
fgout_grid.read_fgout_grids_data_pre511()

# Deterime time t0 of first fgout frame, to initialize particles
frameno0 = fgframes[0]
fgout0 = fgout_grid.read_frame(frameno0)
t0 = fgout0.t

x1fg,x2fg,y1fg,y2fg = fgout0.extent_edges
fgout_extent = [y1fg,y2fg,x2fg,x1fg]  # for rotated
print('fgout0.extent_edges = ', fgout0.extent_edges)
print('fgout_extent = ', fgout_extent)

debris = RigidMotion.DebrisObject()
debris.L = 4 * [0.6]
debris.phi = [0, pi/2, pi/2, pi/2]
debris.height = 0.4
debris.bottom_area = debris.L[0]*debris.L[1]  # assuming rectangle
debris.face_width = debris.L[0]  # assuming square
debris.z0 = [34.34, 0.22, 0]
debris.friction_static = 0.35
debris.friction_kinetic = 0.25
debris.advect = False
mass = 14.5 # kg
debris.rho = mass / (debris.height * debris.bottom_area)
print('Draft = %.2fm' % debris.draft)



t0 = 34.
nsteps = 50 #221
dt = 0.1
corner_paths = RigidMotion.make_corner_paths_accel(debris,h_fcn,u_fcn,v_fcn,
                                                   t0,dt,nsteps,verbose=True)

corner_paths_2 = None  # no comparison

if 1:
    # ===========
    # plotting

    fig,ax = subplots(figsize=(6,7))

    ax.set_xlim(xlimits)
    ax.set_ylim(ylimits)

    #ax.plot([y1b,y1b,y2b,y2b,y1b], [x1b,x2b,x2b,x1b,x1b], 'g')

    imqoi = 'Depth'

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
        
        fgout_h = fgout0.h
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
        
        t = fgout0.t
        t_str = '%.2f seconds' % t
        title_text = title('%s at t = %s' % (imqoi,t_str))
        
    else:
        # speed
        s_units = 'm/s'
        fgout_s = fgout0.s
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

        t = fgout0.t
        t_str = '%.2f seconds' % t
        title_text = title('%s at t = %s' % (imqoi,t_str))

    # plot debris square:
    
    c = {None:'g', 'static':'r', 'kinetic':'orange'}
    
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
            fgout_s = fgout_u_txy[fgframeno,:,:]
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

if __name__ == '__main__':
        
    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=len(corner_paths), 
                                   interval=200, blit=False)

    fname_mp4 = 'corner_paths.mp4'
    fps = 5
    print('Making mp4...')
    animation_tools.make_mp4(anim, fname_mp4, fps)
