from pylab import *
from scipy.interpolate import RegularGridInterpolator
from clawpack.geoclaw import fgout_tools
import glob


#outdir = '../BM1/SWE2d/_output_dtopo1'

def load_fgout(outdir, fgno=1, fgframenos=None):
    format = 'binary32'  # format of fgout grid output

    print('Looking for output in ',outdir)

    output_format = 'binary32'

    # List of frames to use for making debris paths and animation:
    #fgframenos = range(10,121)
    #fgframenos = range(10,30)

    if fgframenos is None:
        fgoutfiles = glob.glob('%s/fgout%s.t*' % (outdir,str(fgno).zfill(4)))
        fgframenos = range(1,len(fgoutfiles)+1)
        print('fgframenos = ',fgframenos)

    # Instantiate object for reading fgout frames:
    fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, output_format)
    try:
        fgout_grid.read_fgout_grids_data_pre511()
    except:
        print('Trying new format...')
        fgout_grid.read_fgout_grids_data()

    frameno0 = fgframenos[0]
    fgout0 = fgout_grid.read_frame(frameno0)
    x = fgout0.x
    y = fgout0.y

    txy_shape = (len(fgframenos), fgout0.X.shape[0], fgout0.X.shape[1])
    fgout_h_txy = empty(txy_shape)
    fgout_u_txy = empty(txy_shape)
    fgout_v_txy = empty(txy_shape)
    fgout_times = empty(len(fgframenos))

    for k,fgframeno in enumerate(fgframenos):
        fgout = fgout_grid.read_frame(fgframeno)
        fgout_times[k] = fgout.t
        fgout_h_txy[k,:,:] = fgout.h
        fgout_u_txy[k,:,:] = fgout.u
        fgout_v_txy[k,:,:] = fgout.v


    print('\nLoaded %i fgout frames at times: \n     %s' \
            % (len(fgframenos), fgout_times))
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

    return fgout_grid, h_fcn, u_fcn, v_fcn, fgout_times, \
                               fgout_h_txy, fgout_u_txy, fgout_v_txy
