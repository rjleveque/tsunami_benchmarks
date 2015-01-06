from pylab import *

d = loadtxt('/Users/rjl/git/tsunami_benchmarks/nthmp_currents_2015/all_data/4_Seaside_OSU_Model_Lab/incident_wave/Wavegage.txt',skiprows=1)
t = d[:,0]
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

savefig('wavemaker_gaussian.png')
