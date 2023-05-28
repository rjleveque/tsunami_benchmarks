
import sys
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend


from pylab import *
import os,sys
from matplotlib import colors
import matplotlib.animation as animation
from clawpack.visclaw import animation_tools, plottools, geoplot
from clawpack.geoclaw import topotools

sys.path.insert(0,'../../common_python')
import fgout_particles as P


if 1:
    from clawpack.geoclaw import fgout_tools
    graphics_dir = os.path.abspath('../../../graphics')
else:
    # local versions for self-contained directory:
    import fgout_tools
    graphics_dir = './'

#topo = topotools.Topography()
#topo.read('../BM3topo.tt3',3)


fgno = 1
fgout_extent = [-1, 7., -4, 0]
#bgimage = [imread(graphics_dir+'fgout02CT.png'), fgout_extent]
bgimage = None
plot_extent = fgout_extent

color = 'r'
linewidth = 3

# List of frames to use for making debris paths and animation:
fgframes = range(1,181)
#fgframes = range(30,35)
#fgframes = range(250,260)

outdir = '_output'
format = 'binary'  # format of fgout grid output

# Create ClawPlotData object used for reading in fgout frames:
#plotdata = P.make_plotdata(fgno, outdir, format)

fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format) 

fgout100 = fgout_grid.read_frame(100)
Bfine = fgout100.B

# Deterime time t0 of first fgout frame, to initialize particles
frameno0 = fgframes[0]
fgout0 = fgout_grid.read_frame(frameno0)
t0 = fgout0.t


def timeformat(t):
    """
    Convert t in seconds to string in format hh:mm:ss
    """
    from numpy import mod
    hours = int(t/3600.)
    tmin = mod(t,3600.)
    min = int(tmin/60.)
    sec = int(mod(tmin,60.))
    timestr = '%s:%s:%s' % (hours,str(min).zfill(2),str(sec).zfill(2))
    return timestr
    
    
# Initialize debris_paths dictionary and set
# initial debris particle locations (x0,y0) at time t0.
# Require a list of dbnos and each array 
#     debris_paths[dbno]
# in the dictionary is a 2d array with a single row [t0, x0, y0] to start.

debris_paths = {}
grounding_depth = {}
drag_factor = {}

dbnos = []
dbnosA = []
dbnosB = []
# set initial velocities to 0 (or may want to interpolate from fgout0 if t0>0?)
u0 = 0.
v0 = 0.

if 0:
    length = 0.05
    xg1 = np.linspace(1.37 - 3.5*length, 1.37+2.5*length, 7)
    yg1 = np.linspace(-2.15-0.05, -2.15+0.05, 3)
    name = 'FD'
if 0:
    length = 0.05
    xg1 = np.linspace(2.94 - 3.5*length, 2.94+2.5*length, 7)
    yg1 = np.linspace(-1.26-0.05, -1.26+0.05, 3)
    name = 'LD'
if 1:
    length = 0.05
    xg1 = np.linspace(0, 5, 10)
    yg1 = np.linspace(-3.5, -1, 10)
    name = 'array'    

#length = 0.25
#xg1 = [1.37 - 0.5*length]
#yg1 = [-2.15 - 0.05, -2.15 + 0.05]

xgg = []
ygg = []
for xg in xg1:
    for yg in yg1:
        xgg.append(xg)
        ygg.append(yg)
        
if 0:
    xg1 = [2.94 - 0.5*length]
    yg1 = [-1.26 - 0.05, -1.26 + 0.05]

    for xg in xg1:
        for yg in yg1:
            xgg.append(xg)
            ygg.append(yg)
        
grounding_depth_common = 0.025
drag_factor_common = None

for k in range(len(xgg)):
    dbno = k
    db = array([[t0, xgg[k], ygg[k], u0, v0]])
    debris_paths[dbno] = db
    dbnos.append(dbno)
    dbnosA.append(dbno)
    grounding_depth[dbno] = grounding_depth_common
    drag_factor[dbno] = drag_factor_common
    
    # add paired particle for long debris:
    dbno = 100+dbno
    db = array([[t0, xgg[k]+length, ygg[k], u0, v0]])
    debris_paths[dbno] = db
    dbnos.append(dbno)
    dbnosB.append(dbno)
    grounding_depth[dbno] = grounding_depth_common
    drag_factor[dbno] = drag_factor_common 


print('Created %i initial debris particles' % len(dbnos))


                     
# Compute debris path for each particle by using all the fgout frames
# in the list fgframes (first frame should be frameno0 used to set t0 above):

debris_paths = P.make_debris_paths_pairs(fgout_grid, fgframes,  
                                         debris_paths, dbnosA, dbnosB,
                                         drag_factor=drag_factor, 
                                         grounding_depth=grounding_depth)


def make_dbAB(t, debris_paths, dbnosA, dbnosB):
    xdAB = []
    ydAB = []
    
    for k in range(len(dbnosA)):
        dbnoA = dbnosA[k]
        dbnoB = dbnosB[k]
        dbA = debris_paths[dbnoA]
        dbB = debris_paths[dbnoB]
        try:
            j = where(abs(dbA[:,0]-t) < 1e-6)[0].max()
        except:
            print('Did not find path for dbno=%i at t = %.3f' % (dbno,t))
            j = -1
        if j > -1:
            xA = dbA[j,1]
            yA = dbA[j,2]
            xB = dbB[j,1]
            yB = dbB[j,2]
            xdAB = xdAB + [xA,xB,nan]
            ydAB = ydAB + [yA,yB,nan]
            #import pdb; pdb.set_trace()
    xdAB = array(xdAB)
    ydAB = array(ydAB)
    return xdAB,ydAB



# First initialize plot with data from initial frame,
# do this in a way that returns an object for each plot attribute that
# will need to be changed in subsequent frames.
# In tis case, the color image of the water depth, plots of particles, and 
# title (which includes the time) will change.
# The background image, colorbar, etc. do not change.

fgout = fgout0
ylat = fgout.Y.mean()  # for aspect ratio of plots

fig,ax = subplots(figsize=(12,7))
if bgimage:
    ax.imshow(bgimage[0],extent=bgimage[1])
    
ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])


s_units = 'm/s'
if s_units == 'knots':
    s = fgout.s * 1.9438  # convert m/s to knots
    bounds_speed = np.array([1e-3,3,4,6,9,12])  # knots
else:
    s = fgout.s
    bounds_speed = np.array([1e-3,0.05,0.1,0.2,0.3,0.4])  # m/s


a = 0.4
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

#contourf(fgout.X, fgout.Y, fgout.B, levels=[0.04, 0.2], colors=[[1,1,.7]])
#Bdry = where(fgout.h < 0.001, fgout.B, nan)
#pcolormesh(fgout.X, fgout.Y, Bdry, cmap=geoplot.land1_colormap)
im = imshow(flipud(s.T), extent=fgout_extent,
            cmap=cmap_speed, norm=norm_speed)
cb = colorbar(im, extend='max', shrink=0.7)
cb.set_label(s_units)
contour(fgout.X, fgout.Y, Bfine, [0], colors='g', linewidths=1)
contour(fgout.X, fgout.Y, Bfine, [0.02,0.04,0.08], colors='g', linewidths=0.5)


ax.set_aspect(1./cos(ylat*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])

t = fgout.t
t_str = timeformat(t)
title_text = title('Water speed at t = %.1f seconds' % t)

xdAB,ydAB = make_dbAB(t, debris_paths, dbnosA, dbnosB)
pairs, = ax.plot(xdAB, ydAB, color=color, linestyle='-', linewidth=linewidth)
print('+++ pairs = ',pairs)

# The function update below should have arguments num (for the frame number)
# plus things listed here in fargs.

fargs = (im,pairs,title_text)

# fargs should be initialized above and are the plot Artist objects 
# whose data change from one frame to the next.


def update(num, im, pairs, title_text):
    
    fgframe = fgframes[num]
    # note: uses fgframes to specify fgout frames to use
    
    # Read fgout data for this frame:
    #fgout = P.read_fgout_frame(fgno, fgframe, plotdata)
    fgout = fgout_grid.read_frame(fgframe)
    
    # Reset the plot objects that need to change from previous frame:

    # title:
    t = fgout.t        
    t_str = timeformat(t)
    title_text.set_text('Water speed at t = %.1f seconds' % t)
    
    # color image:
    im.set_data(flipud(fgout.s.T))

    # particle locations:
    
    xdAB,ydAB = make_dbAB(t, debris_paths, dbnosA, dbnosB)
    pairs.set_data(xdAB, ydAB)        
        
    # must now return all the objects listed in fargs:
    return im,pairs,title_text

print('Making anim...')
anim = animation.FuncAnimation(fig, update,
                               frames=len(fgframes), 
                               fargs=fargs,
                               interval=200, blit=True)

fname_mp4 = '%s_pairs.mp4' % name
fps = 5
print('Making mp4...')
animation_tools.make_mp4(anim, fname_mp4, fps)
