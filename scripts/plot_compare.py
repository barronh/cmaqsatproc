import PseudoNetCDF as pnc
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import numpy as np
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

opts = dict(norm=mc.Normalize())
ropts = dict(norm=mc.BoundaryNorm([-200, -100, -50, -25, -10, 10, 25, 50, 100, 200], 256), cmap='bwr')

fig, axx = plt.subplots(1, 3, figsize=(12, 4))
exec(open(optpath, 'r').read())

satf = pnc.pncopen(satpath, format='ioapi').copy().eval(
    'VCD = {}'.format(satexpr)
).slice(LAY=0).apply(TSTEP='mean')

modf = pnc.pncopen(modpath, format='ioapi').copy().eval(
    'VCD = {}'.format(modexpr)
).apply(TSTEP='mean')

plt.sca(axx[0])
modf.plot('VCD', plot_kw=opts)

plt.sca(axx[1])
satf.plot('VCD', plot_kw=opts)

plt.sca(axx[2])
rf = modf.copy()
rvar = rf.variables['VCD']
rvar[:] = (modf.variables['VCD'][:] / satf.variables['VCD'][:] - 1) * 100
rvar.long_name = 'VCD_NMB'
rvar.units = '%'
rf.plot('VCD', plot_kw=ropts)

plt.savefig(outpath)
