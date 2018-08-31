import matplotlib.pyplot as pl
import pandas as pd
import numpy as np
from scipy.signal import periodogram

from read_DDH_netcdf import *

pl.ion()
pl.close('all')

def abs_velocity(u, v):
    return (u**2 + v**2)**0.5

class Power_spectrum:
    def __init__(self, array, sampling_freq):
        
        # Calculate the power spectrum (units = un**2, with `un` the units of `array`)
        self.k, self.pxx  = periodogram(array, fs=sampling_freq, scaling='spectrum', axis=0)

        # Vertical mean (TO-DO: height weighted)
        self.pxx_m = np.mean(self.pxx, axis=-1)

        # Cumulative
        self.pxx_c   = np.cumsum(self.pxx  [::-1], axis=0)[::-1]
        self.pxx_m_c = np.cumsum(self.pxx_m[::-1]        )[::-1]

        # Frequency in hours:
        self.k[0] = 1e-16       # Mean == 0
        self.f = 1./self.k/3600.

        # Check-check-check....
        var_numpy = np.var(array[:,0])
        var_spec  = np.sum(self.pxx[1:,0])
        print('Variance (k=0) Numpy={}, Spectrum={}'.format(var_numpy, var_spec))


# Period to study
start = datetime.datetime(year=2017, month=1, day=1, hour=0)
end   = datetime.datetime(year=2017, month=9, day=1, hour=0)

# DDH Domain and height:
iloc = 7# +24
kmax = 21   # 21 = <1000 m
kmax = 8    # 8  = <200 m

# Harmonie data from DDH files:
if 'data' not in locals():
    #path  = '/scratch/ms/nl/nkbs/DOWA/LES_forcing/'
    path  = '/nobackup/users/stratum/DOWA/LES_forcing/'

    variables = ['z', 'u', 'v', 'T', 'dtu_dyn', 'dtv_dyn', 'dtT_dyn', 'dtu_phy', 'dtv_phy', 'dtT_phy',]
    data      = read_DDH_netcdf(start, end, path, variables)
    thour     = np.array((data.time - data.time[0]), dtype=np.float) * 1e-9 / 3600.

# Power spectra:
if 'su' not in locals():
    su = Power_spectrum(data['u'][:,iloc,:kmax].values, 1./600.)
    sv = Power_spectrum(data['v'][:,iloc,:kmax].values, 1./600.)
    sT = Power_spectrum(data['T'][:,iloc,:kmax].values, 1./600.)

    sdu = Power_spectrum(data['dtu_dyn'][:,iloc,:kmax].values, 1./600.)
    sdv = Power_spectrum(data['dtv_dyn'][:,iloc,:kmax].values, 1./600.)
    sdT = Power_spectrum(data['dtT_dyn'][:,iloc,:kmax].values, 1./600.)

    spu = Power_spectrum(data['dtu_phy'][:,iloc,:kmax].values, 1./600.)
    spv = Power_spectrum(data['dtv_phy'][:,iloc,:kmax].values, 1./600.)
    spT = Power_spectrum(data['dtT_phy'][:,iloc,:kmax].values, 1./600.)

def time_labels():
    ax = pl.gca()
    yext = ax.get_ybound()

    times  = [1,24,24*30]
    labels = ['hour', 'day', 'month']

    for h in [1,24,24*30]:
        pl.plot([h,h], yext, 'k:')


pl.figure(figsize=(10,6))
pl.subplot(231)
pl.loglog(su.f[1:], su.pxx_m[1:], linewidth=0.5)
time_labels()
pl.xlabel('f (h)')
pl.ylabel('P_uu (m2 s-2)')

pl.subplot(232)
pl.loglog(sv.f[1:], sv.pxx_m[1:], linewidth=0.5)
time_labels()
pl.xlabel('f (h)')
pl.ylabel('P_vv (m2 s-2)')

pl.subplot(233)
pl.loglog(sT.f[1:], sT.pxx_m[1:], linewidth=0.5)
time_labels()
pl.xlabel('f (h)')
pl.ylabel('P_TT (K2)')

pl.subplot(234)
pl.loglog(su.f[1:], su.pxx_m_c[1:]**0.5, linewidth=0.5)
time_labels()
pl.xlabel('f (h)')
pl.ylabel('cumsum(P_uu) (m2 s-2)')

pl.subplot(235)
pl.loglog(sv.f[1:], sv.pxx_m_c[1:]**0.5, linewidth=0.5)
time_labels()
pl.xlabel('f (h)')
pl.ylabel('cumsum(P_vv) (m2 s-2)')

pl.subplot(236)
pl.loglog(sT.f[1:], sT.pxx_m_c[1:]**0.5, linewidth=0.5)
time_labels()
pl.xlabel('f (h)')
pl.ylabel('cumsum(P_TT) (m2 s-2)')

pl.tight_layout()


pl.figure(figsize=(10,6))
pl.subplot(231)
pl.loglog(sdu.f[1:], sdu.pxx_m[1:], linewidth=0.5)
time_labels()

pl.subplot(232)
pl.loglog(sdv.f[1:], sdv.pxx_m[1:], linewidth=0.5)
time_labels()

pl.subplot(233)
pl.loglog(sdT.f[1:], sdT.pxx_m[1:], linewidth=0.5)
time_labels()

pl.subplot(234)
pl.loglog(sdu.f[1:], sdu.pxx_m_c[1:], linewidth=0.5)
time_labels()

pl.subplot(235)
pl.loglog(sdv.f[1:], sdv.pxx_m_c[1:], linewidth=0.5)
time_labels()

pl.subplot(236)
pl.loglog(sdT.f[1:], sdT.pxx_m_c[1:], linewidth=0.5)
time_labels()


