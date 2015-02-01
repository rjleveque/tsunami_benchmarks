
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from clawpack.geoclaw import topotools


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    from clawpack.visclaw import colormaps, geoplot
    from numpy import linspace

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=False)
    
    def fixup(current_data):
        import pylab
        addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Surface at %4.2f hours' % t, fontsize=20)
        pylab.ticklabel_format(format='plain',useOffset=False)
        mean_lat = 19.7
        pylab.gca().set_aspect(1.0 / pylab.cos(pylab.pi / 180.0 * mean_lat))
        #pylab.xticks(fontsize=15)
        #pylab.yticks(fontsize=15)

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    #plotaxes.scaled = True

    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.2
    plotitem.pcolor_cmax = 0.2
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [1,1,0]
    plotitem.patchedges_show = 1
    #plotaxes.xlimits = [202., 206.]
    #plotaxes.ylimits = [19., 21.]

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = linspace(-2000,0,5)
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figure for imshow plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='imshow', figno=1)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('imshow')
    plotaxes.title = 'Surface'
    #plotaxes.scaled = True

    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    # plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = geoplot.tsunami_colormap
    plotitem.imshow_cmin = -0.2
    plotitem.imshow_cmax = 0.2
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = linspace(-2000,0,5)
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figure for zoom plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='zoom', figno=2)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('imshow')
    plotaxes.title = 'Surface'
    #plotaxes.scaled = True

    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    # plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = geoplot.tsunami_colormap
    plotitem.imshow_cmin = -0.2
    plotitem.imshow_cmax = 0.2
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1
    plotaxes.xlimits = [204.8, 205.]
    plotaxes.ylimits = [19.7, 19.9]


    #-----------------------------------------
    # Figure for zoom plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='zoom2', figno=3)
    #plotfigure.kwargs = {'figsize':(6,12)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('imshow')
    plotaxes.title = 'Surface'
    #plotaxes.scaled = True

    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    # plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = geoplot.tsunami_colormap
    plotitem.imshow_cmin = -2.0
    plotitem.imshow_cmax = 2.0
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0
    #plotaxes.xlimits = [204.9003, 204.965]
    #plotaxes.ylimits = [19.71, 19.91]
    plotaxes.xlimits = [204.905,204.965]
    plotaxes.ylimits = [19.71, 19.758]


    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = linspace(-2000,0,5)
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0



    #-----------------------------------------
    # Figure for vorticity plot
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name='vorticity', figno=9)
    #plotfigure.show = False
    #plotfigure.kwargs = {'figsize':(16,6)}


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Vorticity'
    plotaxes.scaled = True
    #plotaxes.xlimits = [204.905,204.965]
    plotaxes.xlimits = [204.905,204.95]
    plotaxes.ylimits = [19.71, 19.758]
    #plotaxes.ylimits = [19.71, 19.75]

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = 4
    plotitem.imshow_cmap = colormaps.blue_white_red
    plotitem.imshow_cmin = -500.
    plotitem.imshow_cmax = 500.
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.6
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    #plotitem.afterpatch = plot_quiver

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = [-0.05, -0.01]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,1]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    def fixup(current_data):
        import pylab
        addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Vorticity at %4.2f hours' % t, fontsize=20)
        pylab.ticklabel_format(format='plain',useOffset=False)
        mean_lat = 19.7
        pylab.gca().set_aspect(1.0 / pylab.cos(pylab.pi / 180.0 * mean_lat))
        pylab.xticks(rotation=20)
    plotaxes.afteraxes = fixup

    #-----------------------------------------
    # Figure for velocity plot
    #-----------------------------------------
    
    plotfigure = plotdata.new_plotfigure(name='velocity', figno=10)
    #plotfigure.show = False

    # Set up for axes for velocity
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.axescmd = 'subplot(212)'
    plotaxes.title = 'Velocity'
    plotaxes.scaled = True

    def speed(current_data):
        from pylab import where,sqrt
        q = current_data.q
        h = q[0,:]
        dry_tol = 0.001
        u = where(h>dry_tol, q[1,:]/h, 0.)
        v = where(h>dry_tol, q[2,:]/h, 0.)
        s = sqrt(u**2 + v**2)
        return s


    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = speed
    #plotitem.imshow_cmap = colormaps.white_red
    #plotitem.imshow_cmap = colormaps.yellow_red_blue
    plotitem.imshow_cmap = \
           colormaps.make_colormap({0:[1,1,1],0.5:[0.5,0.5,1],1:[1,0.3,0.3]})
    plotitem.imshow_cmin = 0.
    plotitem.imshow_cmax = 1.
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.6
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]


    #-----------------------------------------
    # Figure for velocity plot
    #-----------------------------------------
    
    plotfigure = plotdata.new_plotfigure(name='v-velocity', figno=12)
    #plotfigure.show = False

    # Set up for axes for velocity
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.axescmd = 'subplot(212)'
    plotaxes.title = 'Velocity'
    plotaxes.scaled = True

    def uvel(current_data):
        from numpy import where, sqrt
        q = current_data.q
        h = q[0,:]
        dry_tol = 0.001
        u = where(h>dry_tol, q[1,:]/h, 0.)
        return u

    def vvel(current_data):
        from numpy import where, sqrt
        q = current_data.q
        h = q[0,:]
        dry_tol = 0.001
        v = where(h>dry_tol, q[2,:]/h, 0.)
        return v

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = vvel
    #plotitem.imshow_cmap = colormaps.white_red
    #plotitem.imshow_cmap = colormaps.yellow_red_blue
    plotitem.imshow_cmap = colormaps.blue_white_red
    plotitem.imshow_cmin = -1.
    plotitem.imshow_cmax = 1.
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.6
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]




    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'
    plotitem.kwargs = {'linewidth':2}


    def add_zeroline(current_data):
        from pylab import plot, legend, xticks, floor, xlim,ylim
        t = current_data.t
        #legend(('surface','topography'),loc='lower left')
        plot(t, 0*t, 'k')
        n = int(floor(t.max()/1800.)) + 2
        xticks([1800*i for i in range(n)],[str(0.5*i) for i in range(n)])
        #xlim(25000,t.max())
        #ylim(-0.5,0.5)
        #print "+++ gaugeno = ",current_data.gaugeno

    #plotaxes.ylimits = [-0.5, 0.5]
    plotaxes.afteraxes = add_zeroline


    plotfigure = plotdata.new_plotfigure(name='Velocities', figno=301, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.title = 'Velocities'
    #plotaxes.afteraxes = add_zeroline

    # Plot velocity as red curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    def speed(current_data):
        from numpy import where, sqrt
        h = current_data.q[0,:]
        h = where(h>0.01, h, 1.e6)
        u = 100. * current_data.q[1,:] / h
        v = 100. * current_data.q[2,:] / h
        s = sqrt(u**2 + v**2)
        return s
    plotitem.plot_var = speed
    plotitem.plotstyle = 'k-'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def uvel(current_data):
        from numpy import where, sqrt
        h = current_data.q[0,:]
        h = where(h>0.01, h, 1.e6)
        u = 100. * current_data.q[1,:] / h
        return u
    plotitem.plot_var = uvel
    plotitem.plotstyle = 'r-'
    plotitem.kwargs = {'linewidth':2}

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def vvel(current_data):
        from numpy import where, sqrt
        h = current_data.q[0,:]
        h = where(h>0.01, h, 1.e6)
        v = 100. * current_data.q[2,:] / h
        return v
    plotitem.plot_var = vvel
    plotitem.plotstyle = 'g-'
    plotitem.kwargs = {'linewidth':2}

    def add_legend(current_data):
        from pylab import legend
        legend(['u','v'],loc='upper left')
        #legend(['Speed','u','v'],loc='upper left')
        add_zeroline(current_data)
    plotaxes.ylimits = [-50,50]
    plotaxes.afteraxes = add_legend

    #-----------------------------------------
    # Figures for fgmax - max values on fixed grids
    #-----------------------------------------
    otherfigure = plotdata.new_otherfigure(name='max speed',
                    fname='speed.png')
    otherfigure = plotdata.new_otherfigure(name='max elevation', fname='zeta.png')



    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'         # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

