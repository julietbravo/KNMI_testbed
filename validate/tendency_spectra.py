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
        self.var_numpy = np.var(array[:,0])
        self.var_spec  = np.sum(self.pxx[1:,0])
        print('Variance (k=0) Numpy={}, Spectrum={}'.format(self.var_numpy, self.var_spec))

class Power_spectra:
    def __init__(self, data, iloc, kmax):

        freq = 1/600.       # 10-minute data...

        # Absolute velocity and tendencies
        U       = abs_velocity(data['u'      ][:,iloc,:kmax].values, data['v'      ][:,iloc,:kmax].values)
        dtU_dyn = abs_velocity(data['dtu_dyn'][:,iloc,:kmax].values, data['dtv_dyn'][:,iloc,:kmax].values)
        dtU_phy = abs_velocity(data['dtu_phy'][:,iloc,:kmax].values, data['dtv_phy'][:,iloc,:kmax].values)
        dtU     = dtU_dyn + dtU_phy

        self.sU = Power_spectrum(U,                              freq)
        #self.su = Power_spectrum(data['u'][:,iloc,:kmax].values, freq)
        #self.sv = Power_spectrum(data['v'][:,iloc,:kmax].values, freq)
        self.sT = Power_spectrum(data['T'][:,iloc,:kmax].values, freq)

        self.sdU = Power_spectrum(dtU_dyn,                              freq)
        #self.sdu = Power_spectrum(data['dtu_dyn'][:,iloc,:kmax].values, freq)
        #self.sdv = Power_spectrum(data['dtv_dyn'][:,iloc,:kmax].values, freq)
        self.sdT = Power_spectrum(data['dtT_dyn'][:,iloc,:kmax].values, freq)

        self.spU = Power_spectrum(dtU_phy,                              freq)
        #self.spu = Power_spectrum(data['dtu_phy'][:,iloc,:kmax].values, freq)
        #self.spv = Power_spectrum(data['dtv_phy'][:,iloc,:kmax].values, freq)
        self.spT = Power_spectrum(data['dtT_phy'][:,iloc,:kmax].values, freq)

        self.stU = Power_spectrum(dtU,                                                                         freq)
        #self.stu = Power_spectrum(data['dtu_phy'][:,iloc,:kmax].values + data['dtu_dyn'][:,iloc,:kmax].values, freq)
        #self.stv = Power_spectrum(data['dtv_phy'][:,iloc,:kmax].values + data['dtv_dyn'][:,iloc,:kmax].values, freq)
        self.stT = Power_spectrum(data['dtT_phy'][:,iloc,:kmax].values + data['dtT_dyn'][:,iloc,:kmax].values, freq)



# Period to study
start = datetime.datetime(year=2017, month=1,  day=1, hour=0)
end   = datetime.datetime(year=2017, month=11, day=1, hour=0)

# DDH Domain and height:
iloc = 0
kmax = 8    # 8  = <200 m, 21 = < 1000m

# Harmonie data from DDH files:
if 'data' not in locals():
    #path  = '/scratch/ms/nl/nkbs/DOWA/LES_forcing/'
    path  = '/nobackup/users/stratum/DOWA/LES_forcing/'

    variables = ['z', 'u', 'v', 'T', 'dtu_dyn', 'dtv_dyn', 'dtT_dyn', 'dtu_phy', 'dtv_phy', 'dtT_phy',]
    data      = read_DDH_netcdf(start, end, path, variables)
    thour     = np.array((data.time - data.time[0]), dtype=np.float) * 1e-9 / 3600.

# Power spectra:
if 'spo' not in locals():
    spo = Power_spectra(data, iloc,    8)
    s10 = Power_spectra(data, iloc+12, 8)
    s30 = Power_spectra(data, iloc+24, 8)


def time_labels():
    ax = pl.gca()
    yext = ax.get_ybound()

    hours  = [1,3,24,24*30]
    labels = ['1h', '3h', 'day', 'month']

    for hour,label in zip(hours, labels):
        pl.plot([hour, hour], yext, 'k', dashes=[1,1])
        pl.text(hour, yext[1], label, rotation=90, va='top', ha='right')

    ax.set_ylim(yext)

if (True):
    pl.figure(figsize=(10,10))
    pl.subplot(321)
    pl.title('U (m s-1)', loc='left')
    pl.loglog(spo.sU.f[1:], spo.sU.pxx_m[1:], linewidth=0.5, color='C1', label='column')
    pl.loglog(s10.sU.f[1:], s10.sU.pxx_m[1:], linewidth=0.5, color='C2', label='10x10 km')
    pl.loglog(s30.sU.f[1:], s30.sU.pxx_m[1:], linewidth=0.5, color='C3', label='30x30 km')
    pl.ylabel('P (m2 s-2)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')
    time_labels()
    pl.legend()

    pl.subplot(322)
    pl.title('T (K)', loc='left')
    pl.loglog(spo.sT.f[1:], spo.sT.pxx_m[1:], linewidth=0.5, color='C1')
    pl.loglog(s10.sT.f[1:], s10.sT.pxx_m[1:], linewidth=0.5, color='C2')
    pl.loglog(s30.sT.f[1:], s30.sT.pxx_m[1:], linewidth=0.5, color='C3')
    pl.ylabel('P (K2)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')
    time_labels()

    pl.subplot(323)
    pl.loglog(spo.sU.f[2:], spo.sU.pxx_m_c[2:], linewidth=1.5, color='C1')
    pl.loglog(s10.sU.f[2:], s10.sU.pxx_m_c[2:], linewidth=1.5, color='C2')
    pl.loglog(s30.sU.f[2:], s30.sU.pxx_m_c[2:], linewidth=1.5, color='C3')
    pl.ylabel('$c\sum$ P (m2 s-2)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')
    time_labels()

    pl.subplot(324)
    pl.loglog(spo.sT.f[2:], spo.sT.pxx_m_c[2:], linewidth=1.5, color='C1')
    pl.loglog(s10.sT.f[2:], s10.sT.pxx_m_c[2:], linewidth=1.5, color='C2')
    pl.loglog(s30.sT.f[2:], s30.sT.pxx_m_c[2:], linewidth=1.5, color='C3')
    pl.ylabel('$c\sum$ P (K2)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')
    time_labels()

    pl.subplot(325)
    pl.semilogx(spo.sT.f[2:], spo.sU.pxx_m_c[2:]/spo.sU.pxx_m_c[2], linewidth=1.5, color='C1')
    pl.semilogx(s10.sT.f[2:], s10.sU.pxx_m_c[2:]/s10.sU.pxx_m_c[2], linewidth=1.5, color='C2')
    pl.semilogx(s30.sT.f[2:], s30.sU.pxx_m_c[2:]/s30.sU.pxx_m_c[2], linewidth=1.5, color='C3')
    pl.ylabel('$c\sum$ P / total (-)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')
    pl.ylim(0,1)
    time_labels()

    pl.subplot(326)
    pl.semilogx(spo.sT.f[2:], spo.sT.pxx_m_c[2:]/spo.sT.pxx_m_c[2], linewidth=1.5, color='C1')
    pl.semilogx(s10.sT.f[2:], s10.sT.pxx_m_c[2:]/s10.sT.pxx_m_c[2], linewidth=1.5, color='C2')
    pl.semilogx(s30.sT.f[2:], s30.sT.pxx_m_c[2:]/s30.sT.pxx_m_c[2], linewidth=1.5, color='C3')
    pl.ylabel('$c\sum$ P / total (-)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')
    pl.ylim(0,1)
    time_labels()

    pl.tight_layout()

    pl.savefig('power_spectrum_U_T_{}.pdf'.format(iloc))

if (True):

    for size,spec in zip([1,10,30],[spo, s10, s30]):

        pl.figure(figsize=(10,10))
        pl.subplot(321)
        pl.title('$\partial_t U$ ({0:}) (m s-2)'.format(size), loc='left')
        pl.loglog(spec.sdU.f[1:], spec.sdU.pxx_m[1:], linewidth=0.5, color='C1', label='dynamics')
        pl.loglog(spec.spU.f[1:], spec.spU.pxx_m[1:], linewidth=0.5, color='C2', label='physics')
        pl.loglog(spec.stU.f[1:], spec.stU.pxx_m[1:], linewidth=0.5, color='C3', label='total')
        pl.ylabel('P (m2 s-4)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.subplot(322)
        pl.title('$\partial_t T$ (K s-1)', loc='left')
        pl.loglog(spec.sdT.f[1:], spec.sdT.pxx_m[1:], linewidth=0.5, color='C1')
        pl.loglog(spec.spT.f[1:], spec.spT.pxx_m[1:], linewidth=0.5, color='C2')
        pl.loglog(spec.stT.f[1:], spec.stT.pxx_m[1:], linewidth=0.5, color='C3')
        pl.ylabel('P (K2)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.subplot(323)
        pl.loglog(spec.sdU.f[1:], spec.sdU.pxx_m_c[1:], linewidth=1.5, color='C1', label='dynamics')
        pl.loglog(spec.spU.f[1:], spec.spU.pxx_m_c[1:], linewidth=1.0, color='C2', label='physics')
        pl.loglog(spec.stU.f[1:], spec.stU.pxx_m_c[1:], linewidth=1.0, color='C3', label='total')
        pl.ylabel('$c\sum$ P (m2 s-4)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()
        pl.legend(loc=7)

        pl.subplot(324)
        pl.loglog(spec.sdT.f[1:], spec.sdT.pxx_m_c[1:], linewidth=1.5, color='C1')
        pl.loglog(spec.spT.f[1:], spec.spT.pxx_m_c[1:], linewidth=1.0, color='C2')
        pl.loglog(spec.stT.f[1:], spec.stT.pxx_m_c[1:], linewidth=1.0, color='C3')
        pl.ylabel('$c\sum$ P (K2)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.subplot(325)
        pl.semilogx(spec.sdU.f[1:], spec.sdU.pxx_m_c[1:]/spec.sdU.pxx_m_c[1], linewidth=1.5, color='C1')
        pl.semilogx(spec.spU.f[1:], spec.spU.pxx_m_c[1:]/spec.spU.pxx_m_c[1], linewidth=1.0, color='C2')
        pl.semilogx(spec.stU.f[1:], spec.stU.pxx_m_c[1:]/spec.stU.pxx_m_c[1], linewidth=1.0, color='C3')
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0, 1.)
        time_labels()

        pl.subplot(326)
        pl.semilogx(spec.sdT.f[1:], spec.sdT.pxx_m_c[1:]/spec.sdT.pxx_m_c[1], linewidth=1.5, color='C1')
        pl.semilogx(spec.spT.f[1:], spec.spT.pxx_m_c[1:]/spec.spT.pxx_m_c[1], linewidth=1.0, color='C2')
        pl.semilogx(spec.stT.f[1:], spec.stT.pxx_m_c[1:]/spec.stT.pxx_m_c[1], linewidth=1.0, color='C3')
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0,1)
        time_labels()

        pl.tight_layout()

        pl.savefig('power_spectrum_dtU_dtT_{}_{}.pdf'.format(iloc, size))





















#    pl.subplot(234)
#    pl.loglog(sdu.f[1:], sdu.pxx_m_c[1:]/sdu.pxx_m_c[0], color='C1', linewidth=1.5)
#    pl.loglog(spu.f[1:], spu.pxx_m_c[1:]/spu.pxx_m_c[0], color='C2', linewidth=1.5)
#    pl.loglog(stu.f[1:], stu.pxx_m_c[1:]/stu.pxx_m_c[0], color='C3', linewidth=1.5)
#    time_labels()
#    pl.ylim(1e-2, 1.2)
#    pl.ylabel('norm sum(P) (m2 s-4)')
#    pl.xlabel('frequency (h)')
#    pl.grid(axis='y', which='both')
#
#    pl.subplot(235)
#    pl.loglog(sdv.f[1:], sdv.pxx_m_c[1:]/sdv.pxx_m_c[0], color='C1', linewidth=1.5)
#    pl.loglog(spv.f[1:], spv.pxx_m_c[1:]/spv.pxx_m_c[0], color='C2', linewidth=1.5)
#    pl.loglog(stv.f[1:], stv.pxx_m_c[1:]/stv.pxx_m_c[0], color='C3', linewidth=1.5)
#    time_labels()
#    pl.ylim(1e-2, 1.2)
#    pl.ylabel('norm sum(P) (m2 s-4)')
#    pl.xlabel('frequency (h)')
#    pl.grid(axis='y', which='both')
#
#    pl.subplot(236)
#    pl.loglog(sdT.f[1:], sdT.pxx_m_c[1:]/sdT.pxx_m_c[0], color='C1', linewidth=1.5)
#    pl.loglog(spT.f[1:], spT.pxx_m_c[1:]/spT.pxx_m_c[0], color='C2', linewidth=1.5)
#    pl.loglog(stT.f[1:], stT.pxx_m_c[1:]/stT.pxx_m_c[0], color='C3', linewidth=1.5)
#    time_labels()
#    pl.ylim(1e-2, 1.2)
#    pl.ylabel('norm sum(P) (K2 s-2)')
#    pl.xlabel('frequency (h)')
#    pl.grid(axis='y', which='both')



if (False):
    pl.figure(figsize=(10,6))
    pl.subplot(231)
    pl.title('$\partial_t u$ dynamics', loc='left')
    pl.loglog(sdu.f[1:], sdu.pxx_m[1:], linewidth=0.5, color='C1')
    time_labels()
    pl.ylabel('P (m2 s-4)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')

    pl.subplot(232)
    pl.title('$\partial_t v$ dynamics', loc='left')
    pl.loglog(sdv.f[1:], sdv.pxx_m[1:], linewidth=0.5, color='C1')
    time_labels()
    pl.ylabel('P (m2 s-4)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')

    pl.subplot(233)
    pl.title('$\partial_t T$ dynamics', loc='left')
    pl.loglog(sdT.f[1:], sdT.pxx_m[1:], linewidth=0.5, color='C1')
    time_labels()
    pl.ylabel('P (K2)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')

    pl.subplot(234)
    pl.loglog(sdu.f[1:], sdu.pxx_m_c[1:]/sdu.pxx_m_c[0], color='C1', linewidth=1.5)
    pl.loglog(spu.f[1:], spu.pxx_m_c[1:]/spu.pxx_m_c[0], color='C2', linewidth=1.5)
    pl.loglog(stu.f[1:], stu.pxx_m_c[1:]/stu.pxx_m_c[0], color='C3', linewidth=1.5)
    time_labels()
    pl.ylim(1e-2, 1.2)
    pl.ylabel('norm sum(P) (m2 s-4)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')

    pl.subplot(235)
    pl.loglog(sdv.f[1:], sdv.pxx_m_c[1:]/sdv.pxx_m_c[0], color='C1', linewidth=1.5)
    pl.loglog(spv.f[1:], spv.pxx_m_c[1:]/spv.pxx_m_c[0], color='C2', linewidth=1.5)
    pl.loglog(stv.f[1:], stv.pxx_m_c[1:]/stv.pxx_m_c[0], color='C3', linewidth=1.5)
    time_labels()
    pl.ylim(1e-2, 1.2)
    pl.ylabel('norm sum(P) (m2 s-4)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')

    pl.subplot(236)
    pl.loglog(sdT.f[1:], sdT.pxx_m_c[1:]/sdT.pxx_m_c[0], color='C1', linewidth=1.5)
    pl.loglog(spT.f[1:], spT.pxx_m_c[1:]/spT.pxx_m_c[0], color='C2', linewidth=1.5)
    pl.loglog(stT.f[1:], stT.pxx_m_c[1:]/stT.pxx_m_c[0], color='C3', linewidth=1.5)
    time_labels()
    pl.ylim(1e-2, 1.2)
    pl.ylabel('norm sum(P) (K2 s-2)')
    pl.xlabel('frequency (h)')
    pl.grid(axis='y', which='both')

    pl.tight_layout()
