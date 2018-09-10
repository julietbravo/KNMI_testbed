import matplotlib.pyplot as pl
import numpy as np

from read_DDH_netcdf import *

pl.close('all')


def rmse(v1, v2):
    return np.sqrt(((v1 - v2)**2).mean())

def rmean(array, n):
    return np.convolve(array, np.ones((n,))/n, mode='same')

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

        self.time_h = np.array((self.time - self.time[0]), dtype=np.float) * 1e-9
        self.dt   = self.time_h[1:] - self.time_h[:-1]
        self.time_round = np.arange(0, self.time_h.max()+30, 600.)


# Period to study
start = datetime.datetime(year=2017, month=5, day=1, hour=0)
end   = datetime.datetime(year=2017, month=6, day=1, hour=0)

# Don't re-read the data if executed interactively in IPython (with `-i` option)
iloc = 0+24     # 0+24 = 30x30km FINO1
if 'data' not in locals():
    dt = datetime.timedelta(hours=3)
    data = LES_forcings(start-dt, end, iloc)

# Data is read with padding; find indices corresponding to start and end time
t0 = np.abs(data.time - np.datetime64(start)).argmin()
t1 = np.abs(data.time - np.datetime64(end  )).argmin()

nt = t1-t0
nz = data.data.dims['level']
k  = 0

# Arrays with integrated quantities
u = np.zeros(nt)
v = np.zeros(nt)
T = np.zeros(nt)

invalid = np.zeros(nt, dtype=np.bool)
invalid[:] = False

# Initialise values at t=0
u[0] = data.u[t0,k]
v[0] = data.v[t0,k]
T[0] = data.T[t0,k]

n_av = np.arange(1,144,2)

diff_u = np.zeros(n_av.size)
diff_v = np.zeros(n_av.size)
diff_T = np.zeros(n_av.size)

rmse_u = np.zeros(n_av.size)
rmse_v = np.zeros(n_av.size)
rmse_T = np.zeros(n_av.size)

time_min = np.zeros(n_av.size)

for i in range(n_av.size):

    dtu_dyn = rmean(data.dtu_dyn[:,k], n_av[i])
    dtv_dyn = rmean(data.dtv_dyn[:,k], n_av[i])
    dtT_dyn = rmean(data.dtT_dyn[:,k], n_av[i])

    # Time integrate!
    for t in range(t0,t1-1):    # Index in Harmonie data
        tt = t-t0               # Index in offline arrays

        if (data.time_round[t+1]%10800 < 1):
            # Account for data assimilation step in Harmonie
            u[tt+1] = data.u[t+1,k]
            v[tt+1] = data.v[t+1,k]
            T[tt+1] = data.T[t+1,k]

            invalid[tt+1] = True
        else:
            # Simple Euler forward integration...
            u[tt+1] = u[tt] + (dtu_dyn[t+1] + data.dtu_phy[t+1,k]) * data.dt[t]
            v[tt+1] = v[tt] + (dtv_dyn[t+1] + data.dtv_phy[t+1,k]) * data.dt[t]
            T[tt+1] = T[tt] + (dtT_dyn[t+1] + data.dtT_phy[t+1,k]) * data.dt[t]


    # Mask the data assimilation steps...
    u = np.ma.masked_where(invalid, u)
    v = np.ma.masked_where(invalid, v)
    T = np.ma.masked_where(invalid, T)

    # Statistics...
    diff_u[i] = np.mean(u-data.u[t0:t1,k])
    diff_v[i] = np.mean(v-data.v[t0:t1,k])
    diff_T[i] = np.mean(T-data.T[t0:t1,k])

    rmse_u[i] = rmse(data.u[t0:t1,k], u)
    rmse_v[i] = rmse(data.v[t0:t1,k], v)
    rmse_T[i] = rmse(data.T[t0:t1,k], T)

    time_min[i] = (i+1)*10

    print('diff || u: {0:.4f} m/s, v: {1:.4f} m/s, T: {2:.4f} K'.format(diff_u[i], diff_v[i], diff_T[i]))
    print('RMSE || u: {0:.4f} m/s, v: {1:.4f} m/s, T: {2:.4f} K'.format(rmse_u[i], rmse_v[i], rmse_T[i]))


pl.figure()
pl.subplot(211)
pl.plot(data.time, data.dtu_dyn[:,k], label='Dynamics')
pl.plot(data.time, dtu_dyn, label='Physics')
pl.legend()

pl.subplot(212)
pl.plot(data.time, data.u[:,k], label='Harmonie')
pl.plot(data.time[t0:t1], u, label='offline')
pl.legend()


pl.figure()
pl.subplot(121)
pl.plot(time_min/60., rmse_u, label='u')
pl.plot(time_min/60., rmse_v, label='v')
pl.xlabel('Time average (h)')
pl.ylabel('RMSE (m/s)')
pl.legend()
pl.subplot(122)
pl.plot(time_min/60., rmse_T)
pl.xlabel('Time average (h)')
pl.ylabel('RMSE (K)')





#nt   = data.data.dims['time']
#nz   = data.data.dims['level']
#
#nmax = 17
#tscale = np.zeros(nmax)
#rmse_u = np.zeros(nmax)
#rmse_v = np.zeros(nmax)
#rmse_T = np.zeros(nmax)
#
#for n in range(1,nmax+1):
#
#    u = np.zeros((nt, nz))
#    v = np.zeros((nt, nz))
#    T = np.zeros((nt, nz))
#
#    u[0,:] = data.u[0,:]
#    v[0,:] = data.v[0,:]
#    T[0,:] = data.T[0,:]
#
#    dtu_dyn = np.zeros_like(data.dtu_dyn)
#    dtv_dyn = np.zeros_like(data.dtv_dyn)
#    dtT_dyn = np.zeros_like(data.dtT_dyn)
#
#    fac = 0
#    for t in range(nt-n):
#        if t%n == 0:
#            dtu_dyn[t,:] = data.dtu_dyn[t,:]
#            dtv_dyn[t,:] = data.dtv_dyn[t,:]
#            dtT_dyn[t,:] = data.dtT_dyn[t,:]
#
#            t0  = t
#            fac = 0
#        else:
#            fac += 1/n
#            fac0 = 1-fac
#            fac1 = fac
#
#            dtu_dyn[t,:] = fac0 * data.dtu_dyn[t0,:] +  fac1 * data.dtu_dyn[t0+n,:]
#            dtv_dyn[t,:] = fac0 * data.dtv_dyn[t0,:] +  fac1 * data.dtv_dyn[t0+n,:]
#            dtT_dyn[t,:] = fac0 * data.dtT_dyn[t0,:] +  fac1 * data.dtT_dyn[t0+n,:]
#
#    for t in range(nt-n):
#        # Data assimilation step in Harmonie; this is missing
#        # from the dynamica/physics tendencies, so push the
#        # "forecast" back to the data from Harmonie:
#        if (time_round[t+1]%10800 < 1):
#            u[t+1,:] = data.u[t+1,:]
#            v[t+1,:] = data.v[t+1,:]
#            T[t+1,:] = data.T[t+1,:]
#        else:
#            u[t+1,:] = u[t,:] + (dtu_dyn[t+1,:] + data.dtu_phy[t+1,:]) * dt[t]
#            v[t+1,:] = v[t,:] + (dtv_dyn[t+1,:] + data.dtv_phy[t+1,:]) * dt[t]
#            T[t+1,:] = T[t,:] + (dtT_dyn[t+1,:] + data.dtT_phy[t+1,:]) * dt[t]
#
#    tscale[n-1] = n*10.
#    rmse_u[n-1] = rmse(u[::-n,:], data.u[::-n,:])
#    rmse_v[n-1] = rmse(v[::-n,:], data.v[::-n,:])
#    rmse_T[n-1] = rmse(T[::-n,:], data.T[::-n,:])
#
#
#pl.figure()
#pl.plot(tscale, rmse_u, '-x', label='u')
#pl.plot(tscale, rmse_v, '-x', label='v')
#pl.legend()
#pl.ylabel('RMSE (m s-1)')
#pl.xlabel('forcing time (min)')
#
#pl.figure()
#pl.plot(dtu_dyn[:,0], '-x', label='sub-sampled')
#pl.plot(data.dtu_dyn[:,0])
#
#pl.figure()
#pl.subplot(131)
#pl.plot(data.time, data.u[:,0], label='Harmonie')
#pl.plot(data.time, u[:,0], '-x', ms=3, label='Offline integration')
#pl.legend()
#
#pl.subplot(132)
#pl.plot(data.time, data.v[:,0])
#pl.plot(data.time, v[:,0], '-x')
#
#pl.subplot(133)
#pl.plot(data.time, data.T[:,0])
#pl.plot(data.time, T[:,0], '-x')
