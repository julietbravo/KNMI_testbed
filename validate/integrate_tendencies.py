import matplotlib.pyplot as pl
import numpy as np

from read_DDH_netcdf import *

pl.close('all')

# Period to study
start = datetime.datetime(year=2017, month=5, day=1, hour=0)
end   = datetime.datetime(year=2017, month=9, day=1, hour=0)

def rmse(v1, v2):
    return np.sqrt(((v1 - v2)**2).mean())

class LES_forcings:
    def __init__(self, start, end, iloc):
        path      = '/nobackup/users/stratum/DOWA/LES_forcing/'
        variables = ['time', 'z', 'u', 'v', 'T', 'dtu_dyn', 'dtv_dyn', 'dtT_dyn', 'dtu_phy', 'dtv_phy', 'dtT_phy',]
        self.data = read_DDH_netcdf(start, end, path, variables)

        for var in variables:
            if var == 'time':
                setattr(self, var, self.data[var][:].values)
            else:
                setattr(self, var, self.data[var][:,iloc,:].values)

# Don't re-read the data if executed interactively in IPython (with `-i` option)
iloc = 0+24
if 'data' not in locals():
    data = LES_forcings(start, end, iloc)

time = np.array((data.time - data.time[0]), dtype=np.float) * 1e-9
dt   = time[1:] - time[:-1]
time_round = np.arange(0, time.max()+30, 600.)

nt   = data.data.dims['time']
nz   = data.data.dims['level']

nmax = 17
tscale = np.zeros(nmax)
rmse_u = np.zeros(nmax)
rmse_v = np.zeros(nmax)
rmse_T = np.zeros(nmax)

for n in range(1,nmax+1):

    u = np.zeros((nt, nz))
    v = np.zeros((nt, nz))
    T = np.zeros((nt, nz))

    u[0,:] = data.u[0,:]
    v[0,:] = data.v[0,:]
    T[0,:] = data.T[0,:]

    dtu_dyn = np.zeros_like(data.dtu_dyn)
    dtv_dyn = np.zeros_like(data.dtv_dyn)
    dtT_dyn = np.zeros_like(data.dtT_dyn)

    fac = 0
    for t in range(nt-n):
        if t%n == 0:
            dtu_dyn[t,:] = data.dtu_dyn[t,:]
            dtv_dyn[t,:] = data.dtv_dyn[t,:]
            dtT_dyn[t,:] = data.dtT_dyn[t,:]

            t0  = t
            fac = 0
        else:
            fac += 1/n
            fac0 = 1-fac
            fac1 = fac

            dtu_dyn[t,:] = fac0 * data.dtu_dyn[t0,:] +  fac1 * data.dtu_dyn[t0+n,:]
            dtv_dyn[t,:] = fac0 * data.dtv_dyn[t0,:] +  fac1 * data.dtv_dyn[t0+n,:]
            dtT_dyn[t,:] = fac0 * data.dtT_dyn[t0,:] +  fac1 * data.dtT_dyn[t0+n,:]

    for t in range(nt-n):
        # Data assimilation step in Harmonie; this is missing
        # from the dynamica/physics tendencies, so push the
        # "forecast" back to the data from Harmonie:
        if (time_round[t+1]%10800 < 1):
            u[t+1,:] = data.u[t+1,:]
            v[t+1,:] = data.v[t+1,:]
            T[t+1,:] = data.T[t+1,:]
        else:
            u[t+1,:] = u[t,:] + (dtu_dyn[t+1,:] + data.dtu_phy[t+1,:]) * dt[t]
            v[t+1,:] = v[t,:] + (dtv_dyn[t+1,:] + data.dtv_phy[t+1,:]) * dt[t]
            T[t+1,:] = T[t,:] + (dtT_dyn[t+1,:] + data.dtT_phy[t+1,:]) * dt[t]

    tscale[n-1] = n*10.
    rmse_u[n-1] = rmse(u[::-n,:], data.u[::-n,:])
    rmse_v[n-1] = rmse(v[::-n,:], data.v[::-n,:])
    rmse_T[n-1] = rmse(T[::-n,:], data.T[::-n,:])


pl.figure()
pl.plot(tscale, rmse_u, '-x', label='u')
pl.plot(tscale, rmse_v, '-x', label='v')
pl.legend()
pl.ylabel('RMSE (m s-1)')
pl.xlabel('forcing time (min)')

pl.figure()
pl.plot(dtu_dyn[:,0], '-x', label='sub-sampled')
pl.plot(data.dtu_dyn[:,0])

pl.figure()
pl.subplot(131)
pl.plot(data.time, data.u[:,0], label='Harmonie')
pl.plot(data.time, u[:,0], '-x', ms=3, label='Offline integration')
pl.legend()

pl.subplot(132)
pl.plot(data.time, data.v[:,0])
pl.plot(data.time, v[:,0], '-x')

pl.subplot(133)
pl.plot(data.time, data.T[:,0])
pl.plot(data.time, T[:,0], '-x')
