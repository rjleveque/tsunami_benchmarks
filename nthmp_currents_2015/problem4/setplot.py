
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy

datadir = '/Users/rjl/currents_workshop/4_Seaside_OSU_Model_Lab/comparison_data'
try:
    b1 = loadtxt(datadir+'/B1.txt',skiprows=1)
    g201 = reshape(b1,(9000,4))
    b4 = loadtxt(datadir+'/B4.txt',skiprows=1)
    g204 = reshape(b4,(9000,4))
    b6 = loadtxt(datadir+'/B6.txt',skiprows=1)
    g206 = reshape(b6,(9000,4))
    b9 = loadtxt(datadir+'/B9.txt',skiprows=1)
    g209 = reshape(b9,(9000,4))
    data_found = True
except:
    print "warning: error loading comparison data"
    data_found = False



#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'

    def set_drytol(current_data):
        # The drytol parameter is used in masking land and water and
        # affects what color map is used for cells with small water depth h.
        # The cell will be plotted as dry if h < drytol.
        # The best value to use often depends on the application and can
        # be set here (measured in meters):
        current_data.user['drytol'] = 1.e-2

    plotdata.beforeframe = set_drytol

    # To plot gauge locations on imshow or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True, fontsize=8)
    

    #-----------------------------------------
    # Figure for imshow plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='surface', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('surface')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True

    cmap = colormaps.make_colormap({0:[0.5,1,0.5],0.01:[0,1,1], \
                                    0.2:[0,0,1], 0.5:[1,1,0],  0.8:[1,0,0], \
                                    1.:[0.2,0,0]})

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = cmap
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 0.2
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 0.3
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotaxes.xlimits = [0,45]
    plotaxes.ylimits = [-14,14]

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = [0.08]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,1]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    #-----------------------------------------
    # Figure for zoom plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='zoom', figno=2)
    plotfigure.kwargs = {'figsize':(14,10)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('surface')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.afteraxes = addgauges

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = cmap
    plotitem.imshow_cmin = 0.
    plotitem.imshow_cmax = 0.2
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 0.3
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotaxes.xlimits = [32,42]
    plotaxes.ylimits = [-6,3]

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = [0.08]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,1]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    #-----------------------------------------
    # Figure for cross section
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='cross-section', figno=1)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,45]
    plotaxes.ylimits = [-0.2,0.2]
    plotaxes.title = 'Cross section at y=0'

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def xsec(current_data):
        # Return x value and surface eta at this point, along y=0
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q

        ij = find((y <= dy/2.) & (y > -dy/2.))
        x_slice = ravel(x)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        return x_slice, eta_slice


    plotitem.map_2d_to_1d = xsec
    plotitem.plotstyle = 'kx'     ## need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize':3}
    plotitem.amr_show = [1]  # plot on all levels

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def xsec_s(current_data):
        # Return x value and speed at this point, along y=0
        from pylab import find,ravel,where,sqrt
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q
        h = q[0,:]
        dry_tol = 0.001
        u = where(h>dry_tol, q[1,:]/h, 0.)
        v = where(h>dry_tol, q[2,:]/h, 0.)
        s = sqrt(u**2 + v**2)
        s = s / sqrt(9.81/0.97)  # so comparable to eta

        ij = find((y <= dy/2.) & (y > -dy/2.))
        x_slice = ravel(x)[ij]
        s_slice = ravel(s)[ij]
        return x_slice, s_slice

    plotitem.map_2d_to_1d = xsec_s
    plotitem.plotstyle = 'bo'     ## need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize':3}
    plotitem.amr_show = [1]  # plot on all levels

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def xsec_B(current_data):
        # Return x value and B at this point, along y=0
        from pylab import find,ravel,where,sqrt
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        B = eta - h

        ij = find((y <= dy/2.) & (y > -dy/2.))
        x_slice = ravel(x)[ij]
        B_slice = ravel(B)[ij]
        return x_slice, B_slice

    plotitem.map_2d_to_1d = xsec_B
    plotitem.plotstyle = 'g+'     ## need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize':3}
    plotitem.amr_show = [1]  # plot on all levels


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')

    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    #plotaxes.ylimits = [-0.1, 0.2]
    plotaxes.title = 'Surface and speed'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    # Plot speed as red curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def speed(current_data):
        from pylab import where,sqrt
        q = current_data.q
        h = q[0,:]
        dry_tol = 0.001
        u = where(h>dry_tol, q[1,:]/h, 0.)
        v = where(h>dry_tol, q[2,:]/h, 0.)
        s = sqrt(u**2 + v**2)
        #s = s / sqrt(9.81/0.97)  # so comparable to eta
        s = s / 10.  # as in Figure 5 of paper
        return s

    plotitem.plot_var = speed
    plotitem.plotstyle = 'r-'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'


    def add_zeroline(current_data):
        from pylab import plot, legend
        t = current_data.t
        legend(('surface','speed/10'),loc='lower left')
        plot(t, 0*t, 'k')

    plotaxes.afteraxes = add_zeroline


    #-----------------------------------------
    # Figure for patches alone
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='patches', figno=22)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]
    plotaxes.title = 'patches'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.amr_celledges_show = [1,1,0]   
    plotitem.amr_patchedges_show = [1]     

    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
