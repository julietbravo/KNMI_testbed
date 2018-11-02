import matplotlib.pyplot as pl
import xarray as xr
import numpy as np
from itertools import cycle

pl.close('all')

path = '/nobackup/users/stratum/DOWA/nudge_boundary/nudge_boundary_bomex'

class Read_span:
    def __init__(self, path, file, t_start, t_end):
        self.f = xr.open_dataset('{}/{}'.format(path, file))

        # Period to average
        self.t0 = int(np.abs(self.f.time - t_start).argmin())
        self.t1 = int(np.abs(self.f.time - t_end  ).argmin())

        # Time averaged
        try:
            self.um   = self.f['uxz'  ][self.t0:self.t1+1].mean(axis=0)
            self.vm   = self.f['vxz'  ][self.t0:self.t1+1].mean(axis=0)
            self.wm   = self.f['wxz'  ][self.t0:self.t1+1].mean(axis=0)
            self.thlm = self.f['thlxz'][self.t0:self.t1+1].mean(axis=0)
            self.qtm  = self.f['qtxz' ][self.t0:self.t1+1].mean(axis=0)
            self.qlm  = self.f['qlxz' ][self.t0:self.t1+1].mean(axis=0)
        except:
            pass

        self.u2m   = self.f['u2rxz'  ][self.t0:self.t1+1].mean(axis=0)
        self.v2m   = self.f['v2rxz'  ][self.t0:self.t1+1].mean(axis=0)
        self.w2m   = self.f['w2rxz'  ][self.t0:self.t1+1].mean(axis=0)
        self.thl2m = self.f['thl2rxz'][self.t0:self.t1+1].mean(axis=0)
        self.qt2m  = self.f['qt2rxz' ][self.t0:self.t1+1].mean(axis=0)
        self.ql2m  = self.f['ql2rxz' ][self.t0:self.t1+1].mean(axis=0)

        # Space (x) averaged
        self.u2b   = self.u2m  .mean(axis=1)
        self.v2b   = self.v2m  .mean(axis=1)
        self.w2b   = self.w2m  .mean(axis=1)
        self.thl2b = self.thl2m.mean(axis=1)
        self.qt2b  = self.qt2m .mean(axis=1)
        self.ql2b  = self.ql2m .mean(axis=1)

        self.x = self.f['xt']
        self.z = self.f['zt']

ref_sp  = Read_span(path, 'crossxzspan.001.nc', 6*3600, 12*3600)
nud_sp  = Read_span(path, 'crossxzspan.002.nc', 6*3600, 12*3600)
nud_spl = Read_span(path, 'crossxzspan.003.nc', 6*3600, 12*3600)

ref_pr  = xr.open_dataset('{}/profiles.001.nc'.format(path))
nud_pr  = xr.open_dataset('{}/profiles.002.nc'.format(path))
nud_prl = xr.open_dataset('{}/profiles.003.nc'.format(path))


if True:
    # ------------------------------
    # Change variance in streamwise direction
    # ------------------------------
    heights = np.array([100, 350, 850, 1400])
    indexes = np.zeros_like(heights)
    for i,z in enumerate(heights):
        indexes[i] = np.abs(nud_sp.z - (z)).argmin()

    colors = pl.cm.Set1.colors

    x0_ref = 40
    x1_ref = 85

    pl.figure(figsize=(10,6))
    ax=pl.subplot(231)

    for i,k in enumerate(indexes):
        pl.plot(nud_sp.x/1000,  nud_sp.u2m[k,:],  color=colors[i], label=r'$z$ = {0:4.0f} m'.format(heights[i]))
        pl.plot(nud_spl.x/1000, nud_spl.u2m[k,:], color=colors[i], label='', dashes=[4,2])
    for i,k in enumerate(indexes):
        label2 = 'reference' if i == 0 else ''
        pl.plot(x1_ref, ref_sp.u2b[k], '+', color=colors[i], markeredgewidth=1.5, label=label2)
    pl.legend()
    pl.ylim(-0.015, 0.4)
    pl.ylabel(r'$\sigma^2_{u}$ (m$^2$ s$^{-2}$)')

    pl.subplot(232)
    for i,k in enumerate(indexes):
        pl.plot(nud_sp.x/1000,  nud_sp.v2m[k,:], color=colors[i])
        pl.plot(nud_spl.x/1000, nud_spl.v2m[k,:], color=colors[i], dashes=[4,2])
        pl.plot(x1_ref, ref_sp.v2b[k], '+', color=colors[i], markeredgewidth=1.5)
    pl.ylim(-0.015, 0.3)
    pl.ylabel(r'$\sigma^2_{v}$ (m$^2$ s$^{-2}$)')

    pl.subplot(233)
    for i,k in enumerate(indexes):
        pl.plot(nud_sp.x/1000,  nud_sp.w2m[k,:], color=colors[i])
        pl.plot(nud_spl.x/1000, nud_spl.w2m[k,:], color=colors[i], dashes=[4,2])
        pl.plot(x1_ref, ref_sp.w2b[k], '+', color=colors[i], markeredgewidth=1.5)
    pl.ylim(-0.015, 0.4)
    pl.ylabel(r'$\sigma^2_{w}$ (m$^2$ s$^{-2}$)')

    pl.subplot(234)
    for i,k in enumerate(indexes):
        pl.plot(nud_sp.x/1000,  nud_sp.thl2m[k,:], color=colors[i])
        pl.plot(nud_spl.x/1000, nud_spl.thl2m[k,:], color=colors[i], dashes=[4,2])
        pl.plot(x1_ref, ref_sp.thl2b[k], '+', color=colors[i], markeredgewidth=1.5)
    pl.ylim(-0.01, 0.15)
    pl.ylabel(r'$\sigma^2_{\theta}$ (K$^2$)')
    pl.xlabel(r'$x$ (km)')

    pl.subplot(235)
    for i,k in enumerate(indexes):
        pl.plot(nud_sp.x/1000,  nud_sp.qt2m[k,:]*1e6, color=colors[i])
        pl.plot(nud_spl.x/1000, nud_spl.qt2m[k,:]*1e6, color=colors[i], dashes=[4,2])
        pl.plot(x1_ref, ref_sp.qt2b[k]*1e6, '+', color=colors[i], markeredgewidth=1.5)
    pl.ylim(-0.015, 0.5)
    pl.ylabel(r'$\sigma^2_{qt}$ (g$^2$ kg$^{-2}$)')
    pl.xlabel(r'$x$ (km)')

    pl.subplot(236)
    for i,k in enumerate(indexes):
        pl.plot(nud_sp.x/1000, nud_sp.ql2m[k,:]*1e6, color=colors[i])
        pl.plot(nud_spl.x/1000, nud_spl.ql2m[k,:]*1e6, color=colors[i], dashes=[4,2])
        pl.plot(x1_ref, ref_sp.ql2b[k]*1e6, '+', color=colors[i], markeredgewidth=1.5)
    pl.ylim(-0.0002, 0.02)
    pl.ylabel(r'$\sigma^2_{ql}$ (g$^2$ kg$^{-2}$)')
    pl.xlabel(r'$x$ (km)')

    pl.tight_layout()


if False:
    # ------------------------------
    # Vertical profiles variance
    # ------------------------------
    x_locs  = np.array([15, 20, 25, 30, 35, 40, 50])
    indexes = np.zeros_like(x_locs, dtype=np.int)
    for i,x in enumerate(x_locs):
        indexes[i] = np.abs(nud_sp.x - (x*1000)).argmin()

    colors  = pl.cm.tab20.colors   #(np.linspace(0,1,indexes.size))
    lines   = [[1,0], [4,2]]
    lw = 1.5

    pl.figure(figsize=(12,6))
    pl.subplot(131)
    pl.plot(ref_sp.u2b, ref_sp.z, 'k-', label='Reference', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_sp.u2m[:,index], nud_sp.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)], label='x ={0:4.1f} km'.format(x_locs[i]))
    pl.legend()
    pl.xlabel(r'$\sigma^2_{u}$ (m$^2$ s$^{-2}$)')
    pl.ylabel(r'$z$ (m)')

    pl.subplot(132)
    pl.plot(ref_sp.v2b, ref_sp.z, 'k-', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_sp.v2m[:,index], nud_sp.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)])
    pl.xlabel(r'$\sigma^2_{v}$ (m$^2$ s$^{-2}$)')

    pl.subplot(133)
    pl.plot(ref_sp.w2b, ref_sp.z, 'k-', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_sp.w2m[:,index], nud_sp.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)])
    pl.xlabel(r'$\sigma^2_{w}$ (m$^2$ s$^{-2}$)')

    pl.tight_layout()


    pl.figure(figsize=(12,6))
    pl.subplot(131)
    pl.plot(ref_sp.thl2b, ref_sp.z, 'k-', label='Reference', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_sp.thl2m[:,index], nud_sp.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)], label='x ={0:6.2f} km'.format(int(nud_sp.x[index]/1000.)))
    pl.legend()
    pl.xlabel(r'$\sigma^2_{\theta}$ (K$^2$)')
    pl.ylabel(r'$z$ (m)')

    pl.subplot(132)
    pl.plot(ref_sp.qt2b*1e6, ref_sp.z, 'k-', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_sp.qt2m[:,index]*1e6, nud_sp.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)])
    pl.xlabel(r'$\sigma^2_{qt}$ (g$^2$ kg$^{-2}$)')

    pl.subplot(133)
    pl.plot(ref_sp.ql2b*1e6, ref_sp.z, 'k-', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_sp.ql2m[:,index]*1e6, nud_sp.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)])
    pl.xlabel(r'$\sigma^2_{ql}$ (g$^2$ kg$^{-2}$)')
    pl.xlim(0,0.008)

    pl.tight_layout()



if False:
    # ------------------------------
    # Time and spanwise averaged cross-sections
    # ------------------------------

    cmap = pl.cm.terrain_r

    pl.figure()
    pl.subplot(211)
    pl.pcolormesh(nud_sp.x, nud_sp.z, nud_sp.qlm, cmap=cmap)
    pl.colorbar()

    pl.subplot(212)
    pl.pcolormesh(nud_sp.x, nud_sp.z, nud_sp.qtm-nud_sp.qtm.mean(axis=1), cmap=pl.cm.RdBu_r, vmin=-0.0006, vmax=0.0006)
    pl.colorbar()


    pl.figure()
    pl.subplot(311)
    pl.pcolormesh(nud_sp.x, nud_sp.z, nud_sp.u2m, cmap=cmap)
    pl.colorbar()

    pl.subplot(312)
    pl.pcolormesh(nud_sp.x, nud_sp.z, nud_sp.v2m, cmap=cmap)
    pl.colorbar()

    pl.subplot(313)
    pl.pcolormesh(nud_sp.x, nud_sp.z, nud_sp.w2m, cmap=cmap)
    pl.colorbar()

    pl.tight_layout()

    pl.figure()
    pl.subplot(311)
    pl.pcolormesh(nud_sp.x, nud_sp.z, nud_sp.thl2m, cmap=cmap)
    pl.colorbar()

    pl.subplot(312)
    pl.pcolormesh(nud_sp.x, nud_sp.z, nud_sp.qt2m, cmap=cmap)
    pl.colorbar()

    pl.subplot(313)
    pl.pcolormesh(nud_sp.x, nud_sp.z, nud_sp.ql2m, cmap=cmap)
    pl.colorbar()

    pl.tight_layout()

