
from pylab import *
from clawpack.visclaw.data import ClawPlotData



plotdata = ClawPlotData()
plotdata.outdir = '_output'

figure(figsize=(15,4))
for k in range(1,5):
    gnum = 500+k
    subplot(1,4,k)
    g = plotdata.getgauge(gnum)
    plot(g.t, g.q[0,:], 'r')
    u = where(g.q[0,:] > 1e-4, 0.1*g.q[1,:]/g.q[0,:], 0.)
    plot(g.t, u, 'g')
    title('VG-%s' % k)

fname = 'virtual_gauges.png'
savefig(fname)
print "Created ",fname
