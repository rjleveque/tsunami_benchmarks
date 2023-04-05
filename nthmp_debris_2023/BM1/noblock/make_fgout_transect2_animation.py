"""
Make an mp4 animation of fgout grid results. 
This is done in a way that makes the animation quickly and with minimum 
storage required, by making one plot and then defining an update function
that only changes the parts of the plot that change in each frame.
The tuple update_artists contains the list of Artists that must be changed
in update.  Modify this as needed.

"""

from pylab import *
import os
from clawpack.visclaw import plottools, geoplot
from clawpack.visclaw import animation_tools
from matplotlib import animation, colors
from datetime import timedelta
import pickle

from clawpack.geoclaw import fgout_tools
    
with open('../sim_data.pickle','rb') as f:
    sim_data = pickle.load(f)
zeta_fcn = sim_data['zeta_fcn']  # function that interpolates zeta to (t,x)

outdir = '_output'
format = 'binary32'  # format of fgout grid output

if 1:
    fgno = 2  # which fgout grid
    fgframes = range(50,121,2)  # frames of fgout solution to use in animation
    #fgframes = [60,80]
if 0:
    fgno = 1
    fgframes = range(1,91,1)  # frames of fgout solution to use in animation
    fgframes = [10]
    
#fgframes = [60,75]

figsize = (6,4)

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format) 


# Plot one frame of fgout data and define the Artists that will need to
# be updated in subsequent frames:

fgout1 = fgout_grid.read_frame(fgframes[0])

plot_extent = fgout1.extent_edges

fig = figure(figsize=figsize)
ax1 = subplot(2,1,1)

ax1.set_xlim(plot_extent[:2])
ax1.set_ylim(-1, 0.75)
ax1.grid(True)

#jytrans = 0  # index for y transect, lower edge
#jytrans = 28  # index for y transect through obstacle
ytrans = 1.52
jytrans = where(fgout1.y<ytrans)[0].max()
print('jytrans = %i' % jytrans)



#ax1.plot(fgout1.x,fgout1.B[:,jytrans], 'g-')
# piecewise linear bottom:
xB = array([0,10,17.5,32,43.75])
yB = array([0,0,0.5,1,1]) - 0.9017
ax1.plot(xB, yB, 'g-', label='bottom')

eta = ma.masked_where(fgout1.h<0.001, fgout1.eta)

eta_plot, = ax1.plot(fgout1.x,eta[:,jytrans], 'b-', label='surface')

tx = vstack((fgout1.t*ones(fgout1.x.shape), fgout1.x)).T
sim_zeta = zeta_fcn(tx)
sim_zeta_plot, = ax1.plot(fgout1.x, sim_zeta, 'm-', label='sim zeta')

title_text = ax1.set_title('y = %.3f at time %.2f seconds, frame %i' \
                % (ytrans,fgout1.t,fgframes[0]), fontsize=8)

#ax.set_aspect(1./cos(ylat*pi/180.))
ticklabel_format(useOffset=False)
xtick_pos = list(arange(0,50,5))
xtick_labels = ['' for xt in xtick_pos]
#ax1.set_xticklabels([])
ax1.set_xticks(xtick_pos, xtick_labels, fontsize=6)
ax1.set_ylabel('meters', fontsize=8)
ytick_pos = list(arange(-1,0.76,0.25))
ytick_labels = [str(yt) for yt in ytick_pos]
ax1.set_yticks(ytick_pos, ytick_labels, fontsize=6)
ax1.legend(loc='lower right', framealpha=1,fontsize=6)
#xticks(rotation=20)

ax2 = subplot(2,1,2)
uvel = ma.masked_where(fgout1.h<0.001, fgout1.u)
u_plot, = ax2.plot(fgout1.x,uvel[:,jytrans], 'b-', label='velocity')

ax2.legend(loc='lower left', framealpha=1,fontsize=6)
ax2.set_xlim(plot_extent[:2])
ax2.set_ylim(-2,3)
ax2.set_ylabel('m/s', fontsize=8)
ytick_pos = list(arange(-2,3.01,1))
ytick_labels = [str(yt) for yt in ytick_pos]
ax2.set_yticks(ytick_pos, ytick_labels, fontsize=6)
xtick_pos = list(arange(0,50,5))
xtick_labels = [str(xt) for xt in xtick_pos]
ax2.set_xticks(xtick_pos, xtick_labels, fontsize=6)
ax2.grid(True)

# The artists that will be updated for subsequent frames:
update_artists = (eta_plot, u_plot, sim_zeta_plot, title_text)
    
    
        
def update(fgframeno, *update_artists):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    The artists in update_artists must have new data assigned.
    """
    
    fgout = fgout_grid.read_frame(fgframeno)
    print('Updating plot at time %s' % timedelta(seconds=fgout.t))
    
    # unpack update_artists (must agree with definition above):
    eta_plot, u_plot, sim_zeta_plot, title_text = update_artists
        
    # reset title to current time:
    title_text = ax1.set_title('y = %.3f at time %.2f seconds, frame %i' \
                    % (ytrans,fgout.t,fgframeno), fontsize=8)

    # reset surface eta to current state:
    #eta = ma.masked_where(fgout.h<0.001, fgout.eta)
    h = fgout.h[:,jytrans]
    #import pdb; pdb.set_trace()
    jwet = where(h > 0.001)[0].max()
    eta_wet = fgout.eta[:,jytrans]
    eta_wet[(jwet+2):] = nan
    
    eta_plot.set_ydata(eta_wet)
    
    tx = vstack((fgout.t*ones(fgout.x.shape), fgout.x)).T
    sim_zeta = zeta_fcn(tx)
    sim_zeta_plot.set_ydata(sim_zeta)

    u_plot.set_ydata(fgout.u[:,jytrans])
        
    update_artists = (eta_plot, u_plot, sim_zeta_plot, title_text)
    return update_artists

def plot_fgframe(fgframeno):
    """
    Convenience function for plotting one frame.
    But if you use this function in IPython and then try to make the animation,
    it may get into an infinite loop (not sure why).  Close the figure to abort.
    """
    update(fgframeno, *update_artists)
                

def make_anim():
    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=fgframes, 
                                   fargs=update_artists,
                                   interval=200, blit=True)
    return anim

if __name__ == '__main__':

    anim = make_anim()
    
    # Output files:
    name = 'fgout_animation_transect_obstacle'

    fname_mp4 = name + '.mp4'

    fname_html = None
    #fname_html = name + '.html'

    if fname_mp4:
        fps = 5
        print('Making mp4...')
        writer = animation.writers['ffmpeg'](fps=fps)
        anim.save(fname_mp4, writer=writer)
        print("Created %s" % fname_mp4)

    if fname_html:
        # html version:
        animation_tools.make_html(anim, file_name=fname_html, title=name)
    
    
    
    
