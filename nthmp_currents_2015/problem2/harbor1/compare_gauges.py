

from pylab import *
from clawpack.visclaw.data import ClawPlotData
from clawpack.clawutil.data import ClawData

comparison_data_dir = '../../all_data/2_Hilo_2011_Field/comparison_data/'

geodata = ClawData()
geodata.read('geoclaw.data', force=True)
sea_level = geodata.sea_level
print "GeoClaw simulation at sea_level = %g relative to MHW" % sea_level


plotdata = ClawPlotData()

plotdata.outdir = '_output'
#toffset = 7. # hours
toffset = 6.6 # hours fits data at 7760 better

figure(3333, figsize=(12,5))
clf()

# ---  Specified control point ---

CP = loadtxt('../se.dat')
CP_t = CP[:,0] / 60.  # hours since EQ
CP_eta = CP[:,1] + 0.13 + sea_level # adjust to match geoclaw run

g3333 = plotdata.getgauge(3333)
t = g3333.t / 3600. + 7.


plot(CP_t,CP_eta,'k-o',label='Specified')
plot(t, g3333.q[3,:], 'r',label='GeoClaw')
ylabel('meters')
legend(loc='upper right')
xlim(7.0,12.5)
ylim(-1,1.5)
title('Surface elevation at control point')
show()


# ---  Tide Gauge 7760 ---

TG = loadtxt(comparison_data_dir + 'TG_1617760_detided.txt')
TG_t = TG[:,0] / 3600.  # hours since EQ
TG_eta = TG[:,1]

for gaugeno in [7760, 7761, 7762]:
    g= plotdata.getgauge(gaugeno)
    t = g.t / 3600. + toffset
    eta = g.q[3,:] - sea_level   # correct for sea_level to compare with detided

    figure(gaugeno, figsize=(12,5))
    clf()
    plot(TG_t,TG_eta,'k-o',label='Observed')
    plot(t, eta, 'r',label='GeoClaw')
    ylabel('meters')
    legend(loc='upper right')
    xlim(7.0,12.5)
    ylim(-2.5,2)
    title('Tide Gauge gaugeno')
    show()

# ---  ADCP HAI1125 ---

TG = loadtxt(comparison_data_dir + 'HAI1125_detided_harmonic.txt')
TG_t = TG[:,0]  # hours since EQ
TG_u = TG[:,1]
TG_v = TG[:,2]

g = plotdata.getgauge(1125)
t = g.t / 3600. + toffset
u = g.q[1,:] / g.q[0,:]  * 100 # convert to cm/sec
v = g.q[2,:] / g.q[0,:]  * 100 # convert to cm/sec

figure(1125, figsize=(12,10))
clf()
subplot(211)
plot(TG_t,TG_u,'k-o',label='Observed')
plot(t, u, 'r',label='GeoClaw')
ylabel('cm/sec')
legend(loc='upper right')
xlim(7.0,12.5)
ylim(-150,150)
title('HAI1125, E-W velocity')

subplot(212)
plot(TG_t,TG_v,'k-o',label='Observed')
plot(t, v, 'r',label='GeoClaw')
ylabel('cm/sec')
legend(loc='upper right')
xlim(7.0,12.5)
ylim(-150,150)
title('HAI1125, N-S velocity')

show()

# ---  ADCP HAI1126 ---

TG = loadtxt(comparison_data_dir + 'HAI1126_detided_harmonic.txt')
TG_t = TG[:,0]  # hours since EQ
TG_u = TG[:,1]
TG_v = TG[:,2]

g = plotdata.getgauge(1126)
t = g.t / 3600. + toffset
u = g.q[1,:] / g.q[0,:]  * 100 # convert to cm/sec
v = g.q[2,:] / g.q[0,:]  * 100 # convert to cm/sec

figure(1126, figsize=(12,10))
clf()
subplot(211)
plot(TG_t,TG_u,'k-o',label='Observed')
plot(t, u, 'r',label='GeoClaw')
ylabel('cm/sec')
legend(loc='upper right')
xlim(7.0,12.5)
ylim(-150,150)
title('HAI1126, E-W velocity')

subplot(212)
plot(TG_t,TG_v,'k-o',label='Observed')
plot(t, v, 'r',label='GeoClaw')
ylabel('cm/sec')
legend(loc='upper right')
xlim(7.0,12.5)
ylim(-150,150)
title('HAI1126, N-S velocity')

show()

def save_plot(gaugeno):
    fname = 'gauge%s.png' % gaugeno
    figure(gaugeno)
    savefig(fname)
    print "Created ",fname

if 1:
    save_plot(3333)
    save_plot(7760)
    save_plot(7761)
    save_plot(7762)
    save_plot(1125)
    save_plot(1126)
