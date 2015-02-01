"""
Plot fgmax output from GeoClaw runs, assuming points are on a rectangular
grid.

"""

from pylab import *
from numpy import ma
import os
from matplotlib import image

# =========================================================

# Google Earth image for plotting on top of...
# need to change for specific location on coast:

plot_zeta_map = False  # set to false if no image available
#GEmap = image.imread(WAcoast + '/maps/LaPush.png')  
#GEextent = (-124.645,-124.53,47.9,47.946)  # for LaPush
# =========================================================


def make_plots(outdir='_output', plotdir='_plots', gridno=1):

    # Some things that might need to change...

    fgmax_input_file = 'fgmax_grid.txt'
    #fgmax_input_file = 'fgmax%s.txt' % gridno

    plot_zeta = True                     # plot max water depth?
    plot_zeta_times = True               # plot time of max depth?
    plot_arrival_times = True            # plot time of first arrival?
    plot_arrival_times_on_zeta = False   # plot contours of arrival time on depth?

    # The max speed and other quantities may not have been tracked...
    # Requires setting rundata.fgmax_data.num_fgmax_val = 2 or 5 in setrun.py
    plot_speed = True                   # plot max flow speed?
    plot_others = False                  # plot max momentum, mom flux, min depth

    sea_level = 0.  # relative to MHW

    # Contour levels for zeta in filled contour plots:
    #clines_zeta = None  # can set to desired contours of zeta 
    clines_zeta = [0,0.5,1,1.5,2,2.5,3,4]   # meters
    

    # Contour levels for arrival time and time of max zeta:
    #clines_t = None  
    clines_t = linspace(90,240,6)  # minutes

    clines_t_label = clines_t[::2]  # which ones to label 
    clines_t_colors = [.5,.5,.5]    # RGB of color for contour lines
    clines_topo = linspace(0,20,2)  # contours for topography

    # The following quantities might not be available:
    clines_speed = [0,0.5,1,1.5,2,3,4,5,6,8,10]  # contours for speed m/s
    clines_speed = linspace(0,6,13)  # contours for speed m/s
    clines_hs = [0,1,2,3,4,5,6,8,10,20]          # contours for momentum
    clines_hss = [0,2,5,10,25,50,100]            # contours for momentum flux
    clines_min_depth = [0,0.5,1,1.5,2,2.5,3]     # contours for min depth

    # =========================================================

    if not os.path.isdir(outdir):
        raise Exception("Missing directory: %s" % outdir)

    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)

    print "Reading output from ",outdir
    print "Using fgmax input from ",fgmax_input_file

    # read mx and my from the input file:
    try:
        fid = open(fgmax_input_file)
    except:
        raise Exception("cannot open %s" % fgmax_input_file)

    # skip some lines:
    for i in range(6):
        line = fid.readline()

    line = fid.readline().split()
    fid.close()
    mx = int(line[0])
    my = int(line[1])


    fname = outdir + '/fort.FG1.valuemax' 
    print "Reading %s ..." % fname
    try:
        d = loadtxt(fname)
    except:
        raise Exception("*** Cannot read file: %s" % fname)

    ncols = d.shape[1]  
    if ncols not in [6,8,14]:
        print "*** Unexpected number of columns in FG file: ",ncols
    if ncols==6:
        ind_tzeta = 4
        plot_speed = False  # not available
    elif ncols==8:
        ind_s = 4
        ind_tzeta = 5
    else:
        ind_tzeta = 9
        ind_s = 4
        ind_hs = 5
        ind_hss = 6
        ind_minus_h = 7

    x = reshape(d[:,0],(mx,my),order='F')
    y = reshape(d[:,1],(mx,my),order='F')
    y0 = 0.5*(y.min() + y.max())   # mid-latitude for scaling plots

    # Maximum depth:
    h = reshape(d[:,3],(mx,my),order='F')

    # AMR level used for each zeta value:
    level = reshape(d[:,2].astype('int'),(mx,my),order='F')
    
    # Determine topo B at each point from the same level of AMR:
    fname = outdir + '/fort.FG1.aux1' 
    print "Reading %s ..." % fname
    daux = loadtxt(fname)
    topo = []
    nlevels = daux.shape[1]
    for i in range(2,nlevels):
        topoi = reshape(daux[:,i],(mx,my),order='F')
        topoi = ma.masked_where(topoi < -1e50, topoi)
        topo.append(topoi)

    B = ma.masked_where(level==0, topo[0])  # level==0 ==> never updated
    levelmax = level.max()
    for i in range(levelmax):
        B = where(level==i+1, topo[i], B)

    # zeta is defined as maximum depth h onshore and maximum elevation
    # relative to MHW offshore:
    h_or_eta = where(B > 0, h, h+B)
    zeta = h_or_eta

    inundated = logical_and((B>0), (h>0))

    tzeta = reshape(d[:,ind_tzeta],(mx,my),order='F')  # Time maximum h recorded
    tzeta = ma.masked_where(tzeta < -1e50, tzeta)      
    tzeta = ma.masked_where(zeta == 0., tzeta) / 60.  # minutes 

    zeta = ma.masked_where(zeta==0.,zeta)

    if ncols == 8:
        speed = reshape(d[:,ind_s],(mx,my),order='F')
        speed = ma.masked_where(zeta==0.,speed) ###  * 100. # convert to cm/sec

    if ncols == 14:
        speed = reshape(d[:,ind_s],(mx,my),order='F')
        speed = ma.masked_where(zeta==0.,speed) ###  * 100. # convert to cm/sec
        hs = reshape(d[:,ind_hs],(mx,my),order='F')
        hs = ma.masked_where(zeta==0.,hs)
        hss = reshape(d[:,ind_hss],(mx,my),order='F')
        hss = ma.masked_where(zeta==0.,hss)
        minus_h = reshape(d[:,ind_minus_h],(mx,my),order='F')
        minus_h = ma.masked_where(zeta==0.,minus_h)
        min_depth = -minus_h


    # last column is arrival times:
    atimes = reshape(d[:,-1],(mx,my),order='F')
    atimes = ma.masked_where(atimes < -1e50, atimes)  
    atimes = ma.masked_where(zeta == 0., atimes) / 60.  # minutes 

    if clines_zeta is None:
        cmax = zeta.max()
        cmin = zeta.min()
        clines_zeta = linspace(cmin,cmax,10)

    def make_zeta_plot(on_map=False):

        figure(201)
        clf()

        if on_map:
            # plot on map
            imshow(GEmap,extent=GEextent)

        # Plot h or eta along with contours of topo:
        print "max zeta = ", zeta.max()
        colors = discrete_cmap(clines_zeta)
        contourf(x,y,zeta,clines_zeta,colors=colors,extend='max')

        cbar = colorbar()
        cbar.set_ticks(clines_zeta)
        cbar.set_label('meters', fontsize=15)

        # Contours of topo:
        contour(x,y,B,clines_topo,colors='g',linestyles='-')

        if plot_arrival_times_on_zeta:
            # Contours of arrival time
            cs = contour(x,y,atimes,clines_t,colors=clines_t_colors)
            clabel(cs,clines_t_label)

        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        gca().set_aspect(1./cos(y0*pi/180.))  # scale for latitude

        title("Zeta Maximum",fontsize=20)
        if plot_arrival_times_on_zeta:
            title("Zeta Maximum and arrival times",fontsize=15)
        
        if on_map:
            fname = plotdir + '/zeta_map.png' 
        else:
            fname = plotdir + '/zeta.png' 
        savefig(fname)
        print "Created ",fname


    if plot_zeta:
        make_zeta_plot(False)

    if plot_zeta_map:
        make_zeta_plot(True)




    if plot_zeta_times:

        # Plot time max h recorded:
        figure(102)
        clf()

        if clines_t is None:
            clines_t = linspace(tzeta.min(), tzeta.max(), 10)
            

        colors = discrete_cmap_times(clines_t)
        contourf(x,y,tzeta,clines_t,colors=colors,extend='max')
        cbar = colorbar()
        cbar.set_ticks(clines_t)
        cbar.set_label('minutes',fontsize=15)

        # Contours of topo:
        contour(x,y,B,clines_topo,colors='k',linestyles='-')

        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        gca().set_aspect(1./cos(y0*pi/180.))
        title('Time of max zeta', fontsize=20)
        
        fname = plotdir + '/zetatimes.png' 
        savefig(fname)
        print "Created ",fname


    if plot_arrival_times:

        # Plot time max h recorded:
        figure(103)
        clf()
        if clines_t is None:
            clines_t = linspace(atimes.min(), atimes.max(), 10)
        colors = discrete_cmap_times(clines_t)
        contourf(x,y,atimes,clines_t,colors=colors,extend='max')
        cbar = colorbar()
        cbar.set_ticks(clines_t)
        cbar.set_label('minutes',fontsize=15)

        # Contours of topo:
        contour(x,y,B,clines_topo,colors='k',linestyles='-')
        
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        gca().set_aspect(1./cos(y0*pi/180.))

        title('Arrival time', fontsize=20)
        fname = plotdir + '/arrival_times.png' 
        savefig(fname)
        print "Created ",fname


    def plot_variable(name, v, clines, units='m', on_map=False):
        figure()
        clf()

        if on_map:
            # plot on map
            imshow(GEmap,extent=GEextent)

        print "max %s = %s" % (name,v.max())
        colors = discrete_cmap(clines)
        contourf(x,y,v,clines,colors=colors,extend='max')

        cbar = colorbar()
        cbar.set_ticks(clines)
        cbar.set_label(units, fontsize=15)

        # Contours of topo:
        contour(x,y,B,clines_topo,colors='w',linestyles='-')

        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        gca().set_aspect(1./cos(y0*pi/180.))

        title("Maximum %s" % name,fontsize=20)
        
        if on_map:
            fname = plotdir + '/%s_map.png' % name
        else:
            fname = plotdir + '/%s.png' % name
        savefig(fname)
        print "Created ",fname


    if plot_speed:
        for on_map in [False]:
            plot_variable('speed',speed,clines_speed,'m/s',on_map)
    if plot_others:
        for on_map in [False]:
            plot_variable('hs',hs,clines_hs,'m**2 / s',on_map)
            plot_variable('hss',hss,clines_hss, 'm**3 / s**2',on_map)
            plot_variable('min_depth',min_depth,clines_min_depth, 'm', on_map)


def discrete_cmap(clines):
    """
    Construct a discrete color map for the regions between the contour lines
    given in clines. Colors go from turqouise through yellow to red.
    """
    nlines = len(clines)
    n1 = int(floor((nlines-1)/2.))
    n2 = nlines - 1 - n1
    Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
    Red = hstack([linspace(0,0.8,n1), ones(n2)])
    Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
    colors = zip(Red,Green,Blue)
    return colors

def discrete_cmap_times(clines):
    """
    Construct a discrete color map for the regions between the contour lines
    given in clines. For arrival times, colors go from red to turquoise.
    """
    nlines = len(clines)
    n1 = int(floor((nlines-1)/2.))
    n2 = nlines - 1 - n1
    Green = flipud(hstack([linspace(1,1,n1),linspace(1,0,n2)]))
    Red = flipud(hstack([linspace(0,0.8,n1), ones(n2)]))
    Blue = flipud(hstack([linspace(1,0.2,n1), zeros(n2)]))
    colors = zip(Red,Green,Blue)
    return colors

if __name__ == "__main__":
    make_plots()
