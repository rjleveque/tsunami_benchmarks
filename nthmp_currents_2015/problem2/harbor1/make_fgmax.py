"""
Create fgmax_grid.txt and fgmax_transect input files 
"""


from clawpack.geoclaw import fgmax_tools


def make_fgmax_grid():
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2       # will specify a 2d grid of points
    fg.x1 = 204.91 
    fg.x2 = 204.95
    fg.y1 = 19.72
    fg.y2 = 19.75
    fg.dx = 1./3600.    # 1 arcsecond
    fg.tstart_max =  1.5*3600.  # when to start monitoring max values
    fg.tend_max = 1.e10       # when to stop monitoring max values
    fg.dt_check = 10.         # target time (sec) increment between updating 
                               # max values
    fg.min_level_check = 2    # which levels to monitor max on
    fg.arrival_tol = 1.e-2    # tolerance for flagging arrival

    fg.input_file_name = 'fgmax_grid.txt'
    fg.write_input_data()


if __name__ == "__main__":
    make_fgmax_grid()


