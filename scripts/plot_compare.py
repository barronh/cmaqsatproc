import PseudoNetCDF as pnc
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import sys


satpath = sys.argv[1]
satexpr = sys.argv[2]
modpath = sys.argv[3]
modexpr = sys.argv[4]
optpath = sys.argv[5]
outpath = sys.argv[6]
print('sat')
print(satpath)
print(satexpr)
print('mod')
print(modpath)
print(modexpr)
print('opt', optpath)
print('out', outpath)

opts = dict(norm=mc.LogNorm())
ropts = dict(norm=mc.BoundaryNorm(
    [-200, -100, -50, -25, 25, 50, 100, 200],
    256
), cmap='bwr')
copts = dict(orientation='horizontal', pad=0.05)

fig, axx = plt.subplots(
    1, 3, figsize=(12, 4), dpi=200,
    gridspec_kw=dict(left=.05, right=.95, bottom=0.05, top=0.975)
)
plt.setp(axx, facecolor='grey')

satf = pnc.pncopen(satpath, format='ioapi').copy().eval(
    'VCD = {}'.format(satexpr)
).slice(LAY=0).apply(TSTEP='mean')

modf = pnc.pncopen(modpath, format='ioapi').copy().eval(
    'VCD = {}'.format(modexpr)
).apply(TSTEP='mean')

plt.sca(axx[1])
modf.plot('VCD', plot_kw=opts, cbar_kw=copts.copy())

plt.sca(axx[0])
satf.plot('VCD', plot_kw=opts, cbar_kw=copts.copy())

plt.sca(axx[2])
rf = modf.copy()
rvar = rf.variables['VCD']
rvar[:] = (modf.variables['VCD'][:] / satf.variables['VCD'][:] - 1) * 100
rvar.long_name = 'VCD_NMB'
rvar.units = '%'
rf.plot('VCD', plot_kw=ropts, cbar_kw=copts.copy())

exec(open(optpath, 'r').read())

plt.savefig(outpath)
