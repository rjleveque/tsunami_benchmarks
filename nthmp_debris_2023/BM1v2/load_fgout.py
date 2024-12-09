from pylab import *
from scipy.interpolate import RegularGridInterpolator
from clawpack.geoclaw import fgout_tools


outdir = '../BM1/SWE2d/_output_dtopo1'
format = 'binary32'  # format of fgout grid output

print('Looking for output in ',outdir)

output_format = 'binary32'

# List of frames to use for making debris paths and animation:
fgframes = range(10,121)

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(1, outdir, output_format)
fgout_grid.read_fgout_grids_data_pre511()

frameno0 = fgframes[0]
fgout0 = fgout_grid.read_frame(frameno0)
x = fgout0.x
y = fgout0.y

txy_shape = (len(fgframes), fgout0.X.shape[0], fgout0.X.shape[1])
fgout_h_txy = empty(txy_shape)
fgout_u_txy = empty(txy_shape)
fgout_v_txy = empty(txy_shape)
fgout_times = empty(len(fgframes))

for k,fgframeno in enumerate(fgframes):
    fgout = fgout_grid.read_frame(fgframeno)
    fgout_times[k] = fgout.t
    fgout_h_txy[k,:,:] = fgout.h
    fgout_u_txy[k,:,:] = fgout.u
    fgout_v_txy[k,:,:] = fgout.v
    
    
print('\nLoaded %i fgout frames at times: \n     %s' % (len(fgframes), fgout_times))
print('Each frame is on a %i by %i grid with' % (len(x),len(y)))
print('     %g <= x <= %g,  %g <= y <= %g' \
        % (x.min(),x.max(),y.min(),y.max()))

fill_value = 0.

fgout_h_fcn = RegularGridInterpolator((fgout_times,x,y), fgout_h_txy, 
                                       method='linear', bounds_error=False,
                                       fill_value=fill_value)
                                       
def h_fcn(x,y,t):
    val = fgout_h_fcn((array(t), array(x), array(y)))
    try:
        val = float(val)  # if scalar
    except:
        pass
    return val

fgout_u_fcn = RegularGridInterpolator((fgout_times,x,y), fgout_u_txy, 
                                       method='linear', bounds_error=False,
                                       fill_value=fill_value)
                                       
def u_fcn(x,y,t):
    val = fgout_u_fcn((array(t), array(x), array(y)))
    try:
        val = float(val)  # if scalar
    except:
        pass
    return val
    
fgout_v_fcn = RegularGridInterpolator((fgout_times,x,y), fgout_v_txy, 
                                       method='linear', bounds_error=False,
                                       fill_value=fill_value)
                                       
def v_fcn(x,y,t):
    val = fgout_v_fcn((array(t), array(x), array(y)))
    try:
        val = float(val)  # if scalar
    except:
        pass
    return val
