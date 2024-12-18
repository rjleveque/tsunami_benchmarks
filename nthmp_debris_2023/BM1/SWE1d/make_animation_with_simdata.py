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

try:
    from clawpack.geoclaw import geoplot
except:
    print('Could not import geoplot from geoclaw')


from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p
import numpy
    
with open('../sim_data.pickle','rb') as f:
    sim_data = pickle.load(f)
zeta_fcn = sim_data['zeta_fcn']  # function that interpolates zeta to (t,x)
u_vel_fcn = sim_data['u_vel_fcn']  # function that interpolates u_vel to (t,x)

outdir1 = '../SGN1d/_output'
outdir2 = '_output'
format = 'ascii'  # format of fgout grid output

fname_celledges = os.path.abspath('celledges.data')
mapc2p1, mx_edge, xp_edge = make_mapc2p(fname_celledges)


fgframes = range(40,121,1)  # frames of fgout solution to use in animation
#fgframes = [20,40,60]
    

figsize = (6,4)


def read_claw(frameno,outdir):
    fortt_file = '%s/fort.t%s' % (outdir, str(frameno).zfill(4))
    with open(fortt_file,'r') as f:
        line = f.read()
        t1d = float(line.split()[0])
        print('t1d = %.3f' % t1d)
    fortq_file = '%s/fort.q%s' % (outdir, str(frameno).zfill(4))
    q1d = loadtxt(fortq_file, skiprows=5)
    mx1d = q1d.shape[0]
    dx1d = 1./mx1d
    xc1d = arange(dx1d/2, 1, dx1d)
    xp1d = mapc2p1(xc1d)
    h1d = q1d[:,0]
    eta1d = q1d[:,2]
    eta_wet = where(h1d > 0.001, eta1d, nan)
    hu1d = ma.masked_where(h1d < 0.001, q1d[:,1])
    u_wet = hu1d / h1d
    return t1d, xp1d, eta_wet, u_wet

# Plot one frame of fgout data and define the Artists that will need to
# be updated in subsequent frames:

xmin = 0
xmax = 45

fig = figure(figsize=figsize)
ax1 = subplot(2,1,1)

ax1.set_xlim(xmin,xmax)
ax1.set_ylim(-1, 0.75)
ax1.grid(True, lw=0.5)

# piecewise linear bottom:
xB = array([0,10,17.5,32,43.75])
yB = array([0,0,0.5,1,1]) - 0.9017
ax1.plot(xB, yB, 'g-', lw=0.7, label='bottom')

t, x, eta1_wet, u1_wet = read_claw(fgframes[0],outdir1)
eta1_plot, = ax1.plot(x, eta1_wet, 'b-', lw=0.7, label='GeoClaw SGNa')

t, x, eta2_wet, u2_wet = read_claw(fgframes[0],outdir2)
eta2_plot, = ax1.plot(x, eta2_wet, 'c-', lw=0.7, label='GeoClaw SWE')

tx = vstack((t*ones(x.shape), x)).T
sim_zeta = zeta_fcn(tx)
sim_zeta_plot, = ax1.plot(x, sim_zeta, 'r-', lw=0.7, label='sim zeta')

title_text = ax1.set_title('Surface at time %.3f seconds, frame %i' \
                % (t,fgframes[0]), fontsize=8)

#ax.set_aspect(1./cos(ylat*pi/180.))
ticklabel_format(useOffset=False)
xtick_pos = list(arange(0,50,5))
xtick_labels = ['' for xt in xtick_pos]
#ax1.set_xticklabels([])
ax1.set_xticks(xtick_pos, xtick_labels, fontsize=6)
ax1.set_ylabel('meters', fontsize=8)
ytick_pos = list(arange(-1,0.51,0.25))
ytick_labels = [str(yt) for yt in ytick_pos]
ax1.set_yticks(ytick_pos, ytick_labels, fontsize=6)
ax1.legend(loc='lower right', framealpha=1,fontsize=6)
#xticks(rotation=20)

ax2 = subplot(2,1,2)
u1_plot, = ax2.plot(x, u1_wet, 'b-', lw=0.7, label='SGNa velocity')
u2_plot, = ax2.plot(x, u2_wet, 'c-', lw=0.7, label='SWE velocity')

sim_u_vel = u_vel_fcn(tx)
sim_u_plot, = ax2.plot(x, sim_u_vel, 'r-', lw=0.7, label='sim u_vel')

ax2.legend(loc='lower left', framealpha=1,fontsize=6)
ax2.set_xlim(xmin,xmax)
ax2.set_ylim(-2,3)
ax2.set_ylabel('m/s', fontsize=8)
ytick_pos = list(arange(-2,3.01,1))
ytick_labels = [str(yt) for yt in ytick_pos]
ax2.set_yticks(ytick_pos, ytick_labels, fontsize=6)
xtick_pos = list(arange(0,50,5))
xtick_labels = [str(xt) for xt in xtick_pos]
ax2.set_xticks(xtick_pos, xtick_labels, fontsize=6)
ax2.grid(True,lw=0.5)

# The artists that will be updated for subsequent frames:
update_artists = (eta1_plot, u1_plot, eta2_plot, u2_plot, 
                  sim_zeta_plot, sim_u_plot, title_text)
    
#import pdb; pdb.set_trace()
        
def update(fgframeno, *update_artists):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    The artists in update_artists must have new data assigned.
    """
    
    t, x, eta1_wet, u1_wet = read_claw(fgframeno,outdir1)
    t, x, eta2_wet, u2_wet = read_claw(fgframeno,outdir2)

    print('Updating plot at time %s' % timedelta(seconds=t))
    
    # unpack update_artists (must agree with definition above):
    eta1_plot, u1_plot, eta2_plot, u2_plot, \
        sim_zeta_plot, sim_u_plot, title_text = update_artists
        
    # reset title to current time:
    title_text = ax1.set_title('Surface at time %.3f seconds, frame %i' \
                    % (t,fgframeno), fontsize=8)

    # reset surface eta to current state:


    eta1_plot.set_ydata(eta1_wet)
    eta2_plot.set_ydata(eta2_wet)
    
    tx = vstack((t*ones(x.shape), x)).T
    sim_zeta = zeta_fcn(tx)
    sim_zeta_plot.set_ydata(sim_zeta)
    
    sim_u_vel = u_vel_fcn(tx)
    sim_u_plot.set_ydata(sim_u_vel)

    u1_plot.set_ydata(u1_wet)
    u2_plot.set_ydata(u2_wet)
        
    update_artists = (eta1_plot, u1_plot, eta2_plot, u2_plot, 
                      sim_zeta_plot, sim_u_plot, title_text)
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
    name = 'SGN_animation'

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
    
    
    
    
