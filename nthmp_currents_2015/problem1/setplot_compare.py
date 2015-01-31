
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy
from clawpack.visclaw import colormaps
import os

outdir1 = os.path.abspath('_output_manning010_cfl089')
outdir2 = os.path.abspath('_output_manning010_cfl090')


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
    plotdata.outdir = outdir1  # for reading times

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
    # Figure for vorticity plot
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name='vorticity', figno=9)
    #plotfigure.show = False
    #plotfigure.kwargs = {'figsize':(16,6)}
    plotfigure.kwargs = {'figsize':(16,12)}


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.title = 'Vorticity'
    plotaxes.scaled = True
    plotaxes.xlimits = [4.5,9.5]

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.outdir = outdir1
    plotitem.plot_var = 4
    plotitem.imshow_cmap = colormaps.blue_white_red
    plotitem.imshow_cmin = -0.5
    plotitem.imshow_cmax = 0.5
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.6
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    #plotitem.afterpatch = plot_quiver

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.outdir = outdir1
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = [-0.05, -0.01]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,1]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    # Set up for axes for comparison
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.title = 'Vorticity'
    plotaxes.scaled = True
    plotaxes.xlimits = [4.5,9.5]

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.outdir = outdir2
    plotitem.plot_var = 4
    plotitem.imshow_cmap = colormaps.blue_white_red
    plotitem.imshow_cmin = -0.5
    plotitem.imshow_cmax = 0.5
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.6
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    #plotitem.afterpatch = plot_quiver

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.outdir = outdir2
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = [-0.05, -0.01]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,1]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True



    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------

    def add_zeroline(current_data):
        from pylab import plot, legend
        t = current_data.t
        plot(t, 0*t, 'k')

    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')

    plotfigure.kwargs = {'figsize':(14,9)}
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,1)'
    #plotaxes.ylimits = [-0.1, 0.2]
    plotaxes.title = 'Surface'
    plotaxes.afteraxes = add_zeroline

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.outdir = outdir1
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.outdir = outdir2
    plotitem.plot_var = 3
    plotitem.plotstyle = 'r-'
    plotaxes.afteraxes = add_zeroline

    # u-velocity
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,2)'
    #plotaxes.ylimits = [-0.1, 0.2]
    plotaxes.title = 'u-velocity'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def u(current_data):
        from pylab import where,sqrt
        q = current_data.q
        h = q[0,:]
        dry_tol = 0.001
        u = where(h>dry_tol, q[1,:]/h, 0.)
        v = where(h>dry_tol, q[2,:]/h, 0.)
        return u

    plotitem.outdir = outdir1
    plotitem.plot_var = u
    plotitem.plotstyle = 'b-'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.outdir = outdir2
    plotitem.plot_var = u
    plotitem.plotstyle = 'r-'


    # v-velocity
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,3)'
    #plotaxes.ylimits = [-0.1, 0.2]
    plotaxes.title = 'v-velocity'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def v(current_data):
        from pylab import where,sqrt
        q = current_data.q
        h = q[0,:]
        dry_tol = 0.001
        u = where(h>dry_tol, q[1,:]/h, 0.)
        v = where(h>dry_tol, q[2,:]/h, 0.)
        return v

    plotitem.outdir = outdir1
    plotitem.plot_var = v
    plotitem.plotstyle = 'b-'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.outdir = outdir2
    plotitem.plot_var = v
    plotitem.plotstyle = 'r-'


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

    
