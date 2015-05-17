
from pylab import *
from clawpack.visclaw.data import ClawPlotData


datadir = '../all_data/4_Seaside_OSU_Model_Lab/comparison_data'
datadir_other = datadir + '/other'

d = loadtxt(datadir+'/Wavegage.txt',skiprows=1)
t = d[:,0]

def plot_wavemaker():
    dt = 0.02
    wm_disp = 0.47*(d[:,1]-d[0,1])
    wm_speed = (wm_disp[2:] - wm_disp[:-2])/(2.*dt)
    t_speed = t[1:-1]
    
    figure()
    
    plot(t_speed,wm_speed)
    
    beta = 0.25 
    t0 = 14.75
    gaussian = 0.51*exp(-beta*(t_speed - t0)**2)
    plot(t_speed,gaussian,'r')
    xlim(0,40)
    legend(['from data','Gaussian'])
    title('speed of wavemaker')

Wmwg = d[:,2]
wg1 = d[:,3]
wg2 = d[:,4]
wg3 = d[:,5]
wg4 = d[:,6]


# ----

b1 = loadtxt(datadir+'/B1.txt',skiprows=1)
b1a = reshape(b1,(9000,4))
b4 = loadtxt(datadir+'/B4.txt',skiprows=1)
b4a = reshape(b4,(9000,4))
b6 = loadtxt(datadir+'/B6.txt',skiprows=1)
b6a = reshape(b6,(9000,4))
b9 = loadtxt(datadir+'/B9.txt',skiprows=1)
b9a = reshape(b9,(9000,4))



plotdata = ClawPlotData()
#plotdata.outdir = '_output_3'
plotdata.outdir = '_output'


figure(50,figsize=(8,12))
clf()
for gnum,wg in zip([1,2,3,4], [wg1,wg2,wg3,wg4]):
    g = plotdata.getgauge(gnum)
    subplot(4,1,gnum)
    plot(t,wg,'b',label='Measured')
    #plot(g.t, g.q[3,:],'r',label='GeoClaw')
    plot(g.t, g.q[0,:],'r',label='GeoClaw')   # since using wrong gauge_module
    xlim(0,40)
    title('Gauge %s' % gnum)
    ylabel('surface (m)')
legend(loc='upper left')



if 0:
    gauges = {}
    gauges = [0, 1, 2, 3, 4, 5, 101, 102, 103, 104, 105, 106, 107, 108, \
                109, 201, 202, 203, 204, 205, 206, 207, 208, 209]
    for gaugeno in gauges:
        gauges[gaugeno] = plotdata.getgauge(gaugeno)

figure(501); clf()
figure(502); clf()
figure(503); clf()
figure(600); clf()

subp = 0
for bnum,ba in zip([1,4,6,9], [b1a,b4a,b6a,b9a]):
    subp = subp + 1
    g = plotdata.getgauge(200 + bnum)
    figure(200+bnum)
    clf()
    subplot(311)
    plot(ba[:,0],ba[:,1],'b',label='Measured')
    plot(g.t, g.q[1,:],'r',label='GeoClaw')
    xlim(15,40)
    title('Gauge B%s' % bnum)
    ylabel('depth (m)')

    figure(501)
    subplot(4,1,subp)
    plot(ba[:,0],ba[:,1],'b',label='Measured')
    plot(g.t, g.q[1,:],'r',label='GeoClaw')
    xlim(15,40)
    ylabel('depth (m)')
    ymid = array(gca().get_ylim()).mean()
    text(22,ymid,'B%s' % bnum, fontsize=15)
    if subp==1: title('Depth')

    figure(200+bnum)
    subplot(312)
    plot(ba[:,0],ba[:,2],'b')
    #g = plotdata.getgauge(201)
    h = g.q[0,:]
    u = where(h > 0.001, g.q[1,:]/h, 0.)
    v = where(h > 0.001, g.q[2,:]/h, 0.)
    s = sqrt(u**2 + v**2)
    #plot(g.t, s, 'r')
    plot(g.t, u, 'r')
    xlim(15,40)
    #ylabel('speed (m/s)')
    ylabel('u-velocity (m/s)')

    figure(502)
    subplot(4,1,subp)
    plot(ba[:,0],ba[:,2],'b')
    #plot(g.t, s, 'r')
    plot(g.t, u, 'r')
    xlim(15,40)
    ylabel('m/s')
    ymid = array(gca().get_ylim()).mean()
    text(22,ymid,'B%s' % bnum, fontsize=15)
    #if subp==1: title('Speed')
    if subp==1: title('u-velocity')

    figure(200+bnum)
    subplot(313)
    plot(ba[:,0],ba[:,3],'b',label='Measured')
    #g = plotdata.getgauge(201)
    hss = h*s*s
    plot(g.t, hss, 'r',label='GeoClaw')
    xlim(15,40)
    ylabel('mflux (m^3/s^2)')
    legend(loc='upper left')

    figure(503)
    subplot(4,1,subp)
    plot(ba[:,0],ba[:,3],'b',label='Measured')
    plot(g.t, hss, 'r',label='GeoClaw')
    xlim(15,40)
    ylabel('m^3/s^2')
    ymid = array(gca().get_ylim()).mean()
    text(22,ymid,'B%s' % bnum, fontsize=15)
    if subp==1: title('Momentum flux')

    figure(600)
    a = subplot(3,4,subp)
    plot(ba[:,0],ba[:,1],'b',label='Measured')
    plot(g.t, g.q[1,:],'r',label='GeoClaw')
    xlim(20,40)
    ylim(-0.02,0.25)
    a.set_xticklabels([])
    if subp==1: 
        ylabel('depth (m)')
    else:
        a.set_yticklabels([])
    ymid = array(gca().get_ylim()).mean()
    #a.text(17,ymid,'B%s' % bnum, fontsize=15)
    title('Gauge B%s' % bnum, fontsize=15)

    a = subplot(3,4,subp+4)
    plot(ba[:,0],ba[:,2],'b')
    plot(g.t, u, 'r')
    xlim(20,40)
    ylim(-0.05,2.5)
    a.set_xticklabels([])
    if subp==1: 
        ylabel('speed (m/s)')
    else:
        a.set_yticklabels([])
    ymid = array(gca().get_ylim()).mean()
    #a.text(17,ymid,'B%s' % bnum, fontsize=15)

    a = subplot(3,4,subp+8)
    plot(ba[:,0],ba[:,3],'b',label='Measured')
    plot(g.t, hss, 'r',label='GeoClaw')
    xlim(20,40)
    ylim(-0.05,1.0)
    if subp==1: 
        ylabel('mflux (m^3/s^2)')
    else:
        a.set_yticklabels([])
    ymid = array(gca().get_ylim()).mean()
    #a.text(17,ymid,'B%s' % bnum, fontsize=15)


if 0:
    figure(50); fname = 'wg1-4.png'; savefig(fname); print "Saved ",fname
    figure(201); fname = 'B1.png'; savefig(fname); print "Saved ",fname
    figure(204); fname = 'B4.png'; savefig(fname); print "Saved ",fname
    figure(206); fname = 'B6.png'; savefig(fname); print "Saved ",fname
    figure(209); fname = 'B9.png'; savefig(fname); print "Saved ",fname
    figure(501); fname = 'B_depth.png'; savefig(fname); print "Saved ",fname
    figure(502); fname = 'B_velocity.png'; savefig(fname); print "Saved ",fname
    figure(503); fname = 'B_mflux.png'; savefig(fname); print "Saved ",fname
figure(600); fname = 'B_gauges.png'; savefig(fname); print "Saved ",fname

def compare(gaugeno):
    row_num = int(floor(gaugeno/100.))
    #print "row_num = ",row_num
    if row_num == 0:
        return
    elif row_num == 1:
        row = 'A'
    elif row_num == 2:
        row = 'B'
    elif row_num == 3:
        row = 'C'
    elif row_num == 4:
        row = 'D'
    num = str(gaugeno)[2]
    print "Gauge %s%s" % (row,num)
    gauge = "%s%s" % (row,num)
    try:
        gdata = loadtxt(datadir_other + '/Location_%s.txt' % gauge,skiprows=3)
    except:
        gdata = loadtxt(datadir + '/Location_%s.txt' % gauge,skiprows=3)
    gdata = reshape(gdata, (9000,4))

    g = plotdata.getgauge(gaugeno)
    figure(1000)
    clf()
    subplot(311)
    plot(gdata[:,0],gdata[:,1],'b',label='Measured')
    plot(g.t, g.q[1,:],'r',label='GeoClaw')
    xlim(15,40)
    title('Gauge %s' % gauge)
    ylabel('depth (m)')
    
    subplot(312)
    plot(gdata[:,0],gdata[:,2],'b')
    h = g.q[0,:]
    u = where(h > 0.001, g.q[1,:]/h, 0.)
    v = where(h > 0.001, g.q[2,:]/h, 0.)
    s = sqrt(u**2 + v**2)
    plot(g.t, u, 'r')
    xlim(15,40)
    ylabel('u-velocity (m/s)')

    subplot(313)
    plot(gdata[:,0],gdata[:,3],'b',label='Measured')
    hss = h*s*s
    plot(g.t, hss, 'r',label='GeoClaw')
    xlim(15,40)
    ylabel('mflux (m^3/s^2)')
    legend(loc='upper left')

    if 0:
        fname = '%s.png' % gauge
        savefig(fname)
        print "Created ",fname

gaugenos = range(101,110) + range(201,210) + range(301,310) + range(401,405)

for gaugeno in gaugenos:
    compare(gaugeno)
