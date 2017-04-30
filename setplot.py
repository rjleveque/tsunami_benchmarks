""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 
#changed by VL 29042015 to work with GeoClaw v.5.2.2

import numpy as np
import matplotlib.pyplot as plt

from clawpack.geoclaw import topotools

try:
    TG32412 = np.loadtxt('WaveGages.txt')
except:
    print "*** Could not load DART data file"

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot
    from numpy import linspace


    plotdata.clearfigures()  # clear any old figures,axes,items data

    def set_drytol(current_data):
        # The drytol parameter is used in masking land and water and
        # affects what color map is used for cells with small water depth h.
        # The cell will be plotted as dry if h < drytol.
        # The best value to use often depends on the application and can
        # be set here (measured in meters):
#        current_data.user.drytol = 1.e-4
        current_data.user['drytol'] = 1.e-24


    plotdata.beforeframe = set_drytol

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    

    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('imshow')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = [0,5.488]
    plotaxes.ylimits = [0,3.402]


    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = geoplot.tsunami_colormap
    plotitem.imshow_cmin = -0.02
    plotitem.imshow_cmax = 0.02
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1


    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 0.05
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotaxes.xlimits = [0,5.488]
    plotaxes.ylimits = [0,3.402]

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = linspace(-3000,-3000,1)
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface at gauges', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,5.488]
    plotaxes.ylimits = [0,3.402]
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

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
        from pylab import plot, legend, xticks, floor, axis, xlabel
        t = current_data.t 
        gaugeno = current_data.gaugeno

        if gaugeno == 32412:
            try:
                plot(TG32412[:,0], TG32412[:,1], 'r')
                legend(['GeoClaw','Obs'],'lower right')
            except: pass
            axis((0,t.max(),-0.3,0.3))

        plot(t, 0*t, 'k')
        n = int(floor(t.max()/3600.) + 2)
        xticks([3600*i for i in range(n)], ['%i' % i for i in range(n)])
        xlabel('time (hours)')

    plotaxes.afteraxes = add_zeroline


    #-----------------------------------------
    # Figure for imshow plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='imshow', figno=100)
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('imshow')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = [0,5.488]
    plotaxes.ylimits = [0,3.402]

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = geoplot.tsunami_colormap
    plotitem.imshow_cmin = -0.02
    plotitem.imshow_cmax = 0.02
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = [1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 0.05
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotaxes.xlimits = [0,5.488]
    plotaxes.ylimits = [0,3.402]
    

    #-----------------------------------------
    # Figure for imshow plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface and Gauge 1', figno=20)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('imshow')
    plotaxes.show = False
    plotaxes.axescmd = "axes([.1,.5,.8,.4])"
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.afteraxes = addgauges
    plotaxes.xlimits = [0,5.488]
    plotaxes.ylimits = [0,3.402]

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = geoplot.tsunami_colormap
    plotitem.imshow_cmin = -0.03
    plotitem.imshow_cmax = 0.03
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = [1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 0.05
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotaxes.xlimits = [0,5.488]
    plotaxes.ylimits = [0,3.402]


    # Gauge trace:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.show = False
    plotaxes.axescmd = "axes([.1,.1,.8,.3])"
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-0.02, 0.05]
    plotaxes.title = 'Gauge 1'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_gauge_trace')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'
    #plotitem.gaugeno = 1
    #plotitem.gaugenos = 1

    #-----------------------------------------
    # Figure for zoom
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Zoom', figno=10)
    #plotfigure.show = False
    plotfigure.kwargs = {'figsize':[7,7]}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('monai')
    plotaxes.title = 'Monai Valley'
    plotaxes.scaled = True
    #plotaxes.xlimits = [4.0, 5.2]
    #plotaxes.ylimits = [1.3, 2.5]
    plotaxes.xlimits = [4.7, 5.2]
    plotaxes.ylimits = [1.5, 2.2]

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = geoplot.tsunami_colormap
    plotitem.imshow_cmin = -0.02
    plotitem.imshow_cmax = 0.02
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.patchedges_show = [1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 0.05
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = arange(-0.02, 0., .0025)
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True
    
    # Add contour lines of topography:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = arange(0., .2, .01)
    plotitem.amr_contour_colors = ['w']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    # Add dashed contour line for current shoreline
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = 0
    plotitem.contour_levels = [0.002]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'dashed','linewidths':2}
    plotitem.amr_contour_show = [1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True




    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=400, \
                    type='each_gauge')

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,25]
    plotaxes.ylimits = [-0.02, 0.05]
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'
    plotitem.kwargs = {'linewidth':2}

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')

    def gaugetopo(current_data):
        q = current_data.q
        h = q[:,0]
        eta = q[:,3]
        topo = eta - h
        return topo
        
    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'
    def add_zeroline(current_data):
        from pylab import plot, legend
        t = current_data.t
        legend(('surface','topography'),loc='lower left')
        plot(t, 0*t, 'k')
        ##gaugeno = current_data.gaugeno
        #gaugenos = current_data.gaugenos 
        
##        if gaugeno in [5,7,9]:
##            col = (gaugeno-3)/2
#        if gaugenos in [5,7,9]:
#            col = (gaugenos-3)/2
#            plot(labgage[:,0],0.01*labgage[:,col],'r')
#            legend(('GeoClaw','topography','sea level','lab data'),loc='upper left')
#        else:
#            legend(('GeoClaw','topography','sea level'),loc='upper right')
            
        

    plotaxes.afteraxes = add_zeroline


    #-----------------------------------------
    # Figure for grids alone
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='grids', figno=2)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]
    plotaxes.title = 'grids'
    plotaxes.scaled = True

    # Set up for item on these axes:
#    plotitem = plotaxes.new_plotitem(plot_type='2d_grid')
#    plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
#    plotitem.amr_celledges_show = [1,1,0]   
#    plotitem.patchedges_show = [1]     

    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_fignos = 'all'            # list of figures to print    
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?


    return plotdata

    
