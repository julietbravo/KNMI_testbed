import matplotlib.pyplot as pl
import xarray as xr
import pandas as pd
import numpy as np
from scipy.signal import periodogram

pl.ion()
pl.close('all')

# Path to the LES forcings
#path  = '/scratch/ms/nl/nkbs/DOWA/LES_forcing/'
path  = '/nobackup/users/stratum/DOWA/LES_forcing/'

def abs_velocity(u, v):
    return (u**2 + v**2)**0.5


def time_labels():
    ax = pl.gca()
    yext = ax.get_ybound()

    hours  = [1,3,24,24*30,24*365]
    labels = ['1h', '3h', 'day', 'month', 'year']

    for hour,label in zip(hours, labels):
        pl.plot([hour, hour], yext, 'C3', dashes=[1,1])
        pl.text(hour, yext[1], label, rotation=90, va='top', ha='right')

    ax.set_ylim(yext)


def read_harmonie(name):
    print('Reading {}'.format(name))

    # Read the merged (monthly) LES forcings
    files = []
    for y in range(2016,2018):
        for m in range(1,13):
            files.append('{0:}/{1:}_{2:04d}{3:02d}.nc'.format(path,name,y,m))
    ds_hm = xr.open_mfdataset(files)

    return ds_hm


def interpz(data, height, goal_height):
    """
    Fast vectorized interpolation of `data` on time varying `height`
    levels to a single goal_height..
    """

    # 1. Cast to numpy arrays (sometimes xarray can be a b****)
    height = height.values
    data   = data.values

    # Find indices of first model level below/above `goal height`
    k0 = np.absolute(height - goal_height).argmin(axis=1)
    n  = data.shape[0]
    z0 = height[range(n), k0]
    k0[z0 > goal_height] -= 1
    k1 = k0 + 1

    # Calculate interpolation factors
    dzf = height[range(n),k1]-height[range(n),k0]
    dz0 = goal_height-height[range(n),k0]
    dz1 = height[range(n),k1] - goal_height
    f1  = dz0 / dzf
    f0  = dz1 / dzf

    # Interpolate
    return f0*data[range(n),k0] + f1*data[range(n),k1]


class Power_spectrum:
    def __init__(self, array, sampling_freq):

        # Calculate the power spectrum (units = un**2, with `un` the units of `array`)
        self.k, self.pxx  = periodogram(array, fs=sampling_freq, scaling='spectrum', axis=0)

        # Vertical mean (TO-DO: height weighted)
        self.pxx_m = np.mean(self.pxx, axis=-1)

        # Cumulative
        self.pxx_c   = np.cumsum(self.pxx  [::-1], axis=0)[::-1]
        self.pxx_m_c = np.cumsum(self.pxx_m[::-1]        )[::-1]

        ## Frequency in hours:
        self.k[0] = 1e-16       # Mean == 0
        self.f = 1./self.k/3600.

        # Check-check-check....
        self.var_numpy = np.var(array[:,0])
        self.var_spec  = np.sum(self.pxx[1:,0])
        print('Variance (k=0!!) Numpy={0:10.6e}, Spectrum={1:10.6e}'.format(self.var_numpy, self.var_spec))


class Power_spectrum_lev:
    def __init__(self, array, sampling_freq):

        # Calculate the power spectrum (units = un**2, with `un` the units of `array`)
        self.k, self.pxx  = periodogram(array, fs=sampling_freq, scaling='spectrum')

        # Cumulative
        self.pxx_c   = np.cumsum(self.pxx  [::-1])[::-1]

        ## Frequency in hours:
        self.k[0] = 1e-16       # Mean == 0
        self.f = 1./self.k/3600.

        # Check-check-check....
        self.var_numpy = np.var(array)
        self.var_spec  = np.sum(self.pxx)
        print('Variance Numpy={0:10.6e}, Spectrum={1:10.6e}'.format(self.var_numpy, self.var_spec))


class Power_spectra:
    def __init__(self, data, slice):

        freq = 1/600.       # 10-minute data...

        # Absolute velocity and tendencies
        U       = abs_velocity(data['u'      ][:,slice].values, data['v'      ][:,slice].values)
        dtU_dyn = abs_velocity(data['dtu_dyn'][:,slice].values, data['dtv_dyn'][:,slice].values)
        dtU_phy = abs_velocity(data['dtu_phy'][:,slice].values, data['dtv_phy'][:,slice].values)
        dtU     = dtU_dyn + dtU_phy

        self.sU  = Power_spectrum(U, freq)
        self.sT  = Power_spectrum(data['T'][:,slice].values, freq)

        self.sdU = Power_spectrum(dtU_dyn, freq)
        self.sdT = Power_spectrum(data['dtT_dyn'][:,slice].values, freq)

        self.spU = Power_spectrum(dtU_phy, freq)
        self.spT = Power_spectrum(data['dtT_phy'][:,slice].values, freq)

        self.stU = Power_spectrum(dtU, freq)
        self.stT = Power_spectrum(data['dtT_phy'][:,slice].values + data['dtT_dyn'][:,slice].values, freq)


if (True):
    # ----------------------------
    # Spectra: Harmonie
    # ----------------------------

    # Read the LES forcings
    if 'cb_30km' not in locals():
        cb_30km = read_harmonie('Cabauw_30km')
        lb_30km = read_harmonie('Loobos_30km')
        k1_30km = read_harmonie('K13_30km')

    # Create power spectra
    if 'ps_cb' not in locals():
        slice = np.s_[0:8]       # 8=<200m, 21=<1000m

        ps_cb = Power_spectra(cb_30km, slice)
        ps_lb = Power_spectra(lb_30km, slice)
        ps_k1 = Power_spectra(k1_30km, slice)

    if True:
        pl.figure(figsize=(10,10))
        pl.subplot(321)
        pl.title('U (m s-1)', loc='left')
        pl.loglog(ps_cb.sU.f[1:], ps_cb.sU.pxx_m[1:], linewidth=0.5)
        pl.ylabel('P (m2 s-2)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()
        pl.legend()

        pl.subplot(322)
        pl.title('T (K)', loc='left')
        pl.loglog(ps_cb.sT.f[1:], ps_cb.sT.pxx_m[1:], linewidth=0.5)
        pl.ylabel('P (K2)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.subplot(323)
        pl.loglog(ps_cb.sU.f[2:], ps_cb.sU.pxx_m_c[2:], linewidth=1.5, color='k')
        pl.ylabel('$c\sum$ P (m2 s-2)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.subplot(324)
        pl.loglog(ps_cb.sT.f[2:], ps_cb.sT.pxx_m_c[2:], linewidth=1.5, color='k')
        pl.ylabel('$c\sum$ P (K2)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.subplot(325)
        pl.semilogx(ps_cb.sT.f[2:], ps_cb.sU.pxx_m_c[2:]/ps_cb.sU.pxx_m_c[2], linewidth=1.5, color='k')
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0,1)
        time_labels()

        pl.subplot(326)
        pl.semilogx(ps_cb.sT.f[2:], ps_cb.sT.pxx_m_c[2:]/ps_cb.sT.pxx_m_c[2], linewidth=1.5, color='k')
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0,1)
        time_labels()

        pl.tight_layout()

        pl.savefig('power_spectrum_U_T.pdf')


    if True:
        pl.figure(figsize=(10,10))
        pl.subplot(321)
        pl.title('$\partial_t U$ (m s-2)', loc='left')
        pl.loglog(ps_cb.sdU.f[1:], ps_cb.sdU.pxx_m[1:], linewidth=0.5, color='k',  label='dynamics')
        pl.loglog(ps_cb.spU.f[1:], ps_cb.spU.pxx_m[1:], linewidth=0.5, color='C2', label='physics', dashes=[2,2])
        pl.loglog(ps_cb.stU.f[1:], ps_cb.stU.pxx_m[1:], linewidth=0.5, color='C1', label='total', dashes=[4,1])
        pl.ylabel('P (m2 s-4)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.subplot(322)
        pl.title('$\partial_t T$ (K s-1)', loc='left')
        pl.loglog(ps_cb.sdT.f[1:], ps_cb.sdT.pxx_m[1:], linewidth=0.5, color='k')
        pl.loglog(ps_cb.spT.f[1:], ps_cb.spT.pxx_m[1:], linewidth=0.5, color='C2', dashes=[2,2])
        pl.loglog(ps_cb.stT.f[1:], ps_cb.stT.pxx_m[1:], linewidth=0.5, color='C1', dashes=[4,1])
        pl.ylabel('P (K2)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.subplot(323)
        pl.loglog(ps_cb.sdU.f[1:], ps_cb.sdU.pxx_m_c[1:], linewidth=1.5, color='k', label='dynamics')
        pl.loglog(ps_cb.spU.f[1:], ps_cb.spU.pxx_m_c[1:], linewidth=1.5, color='C2', label='physics', dashes=[2,2])
        pl.loglog(ps_cb.stU.f[1:], ps_cb.stU.pxx_m_c[1:], linewidth=1.5, color='C1', label='total', dashes=[4,1])
        pl.ylabel('$c\sum$ P (m2 s-4)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()
        pl.legend(loc=7)

        pl.subplot(324)
        pl.loglog(ps_cb.sdT.f[1:], ps_cb.sdT.pxx_m_c[1:], linewidth=1.5, color='k')
        pl.loglog(ps_cb.spT.f[1:], ps_cb.spT.pxx_m_c[1:], linewidth=1.5, color='C2', dashes=[2,2])
        pl.loglog(ps_cb.stT.f[1:], ps_cb.stT.pxx_m_c[1:], linewidth=1.5, color='C1', dashes=[4,1])
        pl.ylabel('$c\sum$ P (K2)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.subplot(325)
        pl.semilogx(ps_cb.sdU.f[1:], ps_cb.sdU.pxx_m_c[1:]/ps_cb.sdU.pxx_m_c[1], linewidth=1.5, color='k')
        pl.semilogx(ps_cb.spU.f[1:], ps_cb.spU.pxx_m_c[1:]/ps_cb.spU.pxx_m_c[1], linewidth=1.5, color='C2', dashes=[2,2])
        pl.semilogx(ps_cb.stU.f[1:], ps_cb.stU.pxx_m_c[1:]/ps_cb.stU.pxx_m_c[1], linewidth=1.5, color='C1', dashes=[4,1])
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0, 1.)
        time_labels()

        pl.subplot(326)
        pl.semilogx(ps_cb.sdT.f[1:], ps_cb.sdT.pxx_m_c[1:]/ps_cb.sdT.pxx_m_c[1], linewidth=1.5, color='k')
        pl.semilogx(ps_cb.spT.f[1:], ps_cb.spT.pxx_m_c[1:]/ps_cb.spT.pxx_m_c[1], linewidth=1.5, color='C2', dashes=[2,2])
        pl.semilogx(ps_cb.stT.f[1:], ps_cb.stT.pxx_m_c[1:]/ps_cb.stT.pxx_m_c[1], linewidth=1.5, color='C1', dashes=[4,1])
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0,1)
        time_labels()

        pl.tight_layout()

        pl.savefig('power_spectrum_dtU_dtT.pdf')


    if False:
        pl.figure(figsize=(10,10))

        pl.subplot(321)
        pl.title('dynamics, U', loc='left')
        pl.semilogx(ps_cb.sdU.f[1:], ps_cb.sdU.pxx_m_c[1:]/ps_cb.sdU.pxx_m_c[1], linewidth=1, color='k', label='Cabauw')
        pl.semilogx(ps_lb.sdU.f[1:], ps_lb.sdU.pxx_m_c[1:]/ps_lb.sdU.pxx_m_c[1], linewidth=1, color='C2', dashes=[2,2], label='Loobos')
        pl.semilogx(ps_k1.sdU.f[1:], ps_k1.sdU.pxx_m_c[1:]/ps_k1.sdU.pxx_m_c[1], linewidth=1, color='C1', dashes=[4,1], label='K13')
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0, 1.)
        time_labels()
        pl.legend()

        pl.subplot(322)
        pl.title('dynamics, T', loc='left')
        pl.semilogx(ps_cb.sdT.f[1:], ps_cb.sdT.pxx_m_c[1:]/ps_cb.sdT.pxx_m_c[1], linewidth=1, color='k')
        pl.semilogx(ps_lb.sdT.f[1:], ps_lb.sdT.pxx_m_c[1:]/ps_lb.sdT.pxx_m_c[1], linewidth=1, color='C2', dashes=[2,2])
        pl.semilogx(ps_k1.sdT.f[1:], ps_k1.sdT.pxx_m_c[1:]/ps_k1.sdT.pxx_m_c[1], linewidth=1, color='C1', dashes=[4,1])
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0,1)
        time_labels()

        pl.subplot(323)
        pl.title('physics, U', loc='left')
        pl.semilogx(ps_cb.spU.f[1:], ps_cb.spU.pxx_m_c[1:]/ps_cb.spU.pxx_m_c[1], linewidth=1, color='k')
        pl.semilogx(ps_lb.spU.f[1:], ps_lb.spU.pxx_m_c[1:]/ps_lb.spU.pxx_m_c[1], linewidth=1, color='C2', dashes=[2,2])
        pl.semilogx(ps_k1.spU.f[1:], ps_k1.spU.pxx_m_c[1:]/ps_k1.spU.pxx_m_c[1], linewidth=1, color='C1', dashes=[4,1])
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0, 1.)
        time_labels()

        pl.subplot(324)
        pl.title('physics, T', loc='left')
        pl.semilogx(ps_cb.spT.f[1:], ps_cb.spT.pxx_m_c[1:]/ps_cb.spT.pxx_m_c[1], linewidth=1, color='k')
        pl.semilogx(ps_lb.spT.f[1:], ps_lb.spT.pxx_m_c[1:]/ps_lb.spT.pxx_m_c[1], linewidth=1, color='C2', dashes=[2,2])
        pl.semilogx(ps_k1.spT.f[1:], ps_k1.spT.pxx_m_c[1:]/ps_k1.spT.pxx_m_c[1], linewidth=1, color='C1', dashes=[4,1])
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0,1)
        time_labels()

        pl.subplot(325)
        pl.title('total, U', loc='left')
        pl.semilogx(ps_cb.stU.f[1:], ps_cb.stU.pxx_m_c[1:]/ps_cb.stU.pxx_m_c[1], linewidth=1, color='k')
        pl.semilogx(ps_lb.stU.f[1:], ps_lb.stU.pxx_m_c[1:]/ps_lb.stU.pxx_m_c[1], linewidth=1, color='C2', dashes=[2,2])
        pl.semilogx(ps_k1.stU.f[1:], ps_k1.stU.pxx_m_c[1:]/ps_k1.stU.pxx_m_c[1], linewidth=1, color='C1', dashes=[4,1])
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0, 1.)
        time_labels()

        pl.subplot(326)
        pl.title('total, T', loc='left')
        pl.semilogx(ps_cb.stT.f[1:], ps_cb.stT.pxx_m_c[1:]/ps_cb.stT.pxx_m_c[1], linewidth=1, color='k')
        pl.semilogx(ps_lb.stT.f[1:], ps_lb.stT.pxx_m_c[1:]/ps_lb.stT.pxx_m_c[1], linewidth=1, color='C2', dashes=[2,2])
        pl.semilogx(ps_k1.stT.f[1:], ps_k1.stT.pxx_m_c[1:]/ps_k1.stT.pxx_m_c[1], linewidth=1, color='C1', dashes=[4,1])
        pl.ylabel('$c\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0,1)
        time_labels()

        pl.tight_layout()




    # ------------- Supporting figures -------------
    if False:

        levs = np.array([4,14,20])
        z    = np.array([100,500,1000])

        # Montly averaged dynamic temperature tendency
        dTm  = np.zeros((levs.size, 12))
        for k,lev in enumerate(levs):
            df = cb_30km['dtT_dyn'].sel(level=lev).to_dataframe()
            dTm[k,:] = np.squeeze( df.groupby(df.index.month).mean() )

        # Daily cycle dynamic temperature tendency
        dTd_mjja = np.zeros((levs.size, 24))
        dTd_ndjf = np.zeros((levs.size, 24))
        for k,lev in enumerate(levs):
            df   = cb_30km['dtT_dyn'].sel(level=lev).to_dataframe()

            df_s = df.loc[(df.index.month >= 5)&(df.index.month <=8)]
            df_w = df.loc[(df.index.month >= 11)|(df.index.month <=2)]

            dTd_mjja[k,:] = np.squeeze( df_s.groupby(df_s.index.hour).mean() )
            dTd_ndjf[k,:] = np.squeeze( df_w.groupby(df_w.index.hour).mean() )


        pl.figure()
        ax=pl.subplot(121)
        for k in range(levs.size):
            pl.plot(np.arange(1,13), dTm[k,:]*3600, label='z={} m'.format(z[k]))
        pl.legend()
        pl.ylabel('dtT_dyn (K h-1)')
        pl.xlabel('month (-)')
        pl.xlim(1,12)
        ax.set_xticks(np.arange(1,13))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])

        pl.subplot(122)
        pl.plot(np.arange(0,24), dTd_mjja[0,:]*3600, label='z=100 m, MJJA')
        pl.plot(np.arange(0,24), dTd_ndjf[0,:]*3600, label='z=100 m, NDJF')
        pl.legend()
        pl.ylabel('dtT_dyn (K h-1)')
        pl.xlabel('time UTC (h)')
        pl.xlim(0,23)




if (False):
    # ----------------------------
    # Spectra: Harmonie vs Cabauw
    # ----------------------------

    # Read the merged (monthly) LES forcings
    files = []
    for y in range(y0,y1+1):
        for m in range(m0,m1+1):
            files.append('{0:}/Cabauw_column_{1:04d}{2:02d}.nc'.format(path,y,m))
    ds_hm = xr.open_mfdataset(files)

    # Interpolate the model levels to the observations heights of Cabauw
    hm_u200 = interpz(ds_hm['u'], ds_hm['z'], 200)
    hm_v200 = interpz(ds_hm['v'], ds_hm['z'], 200)
    hm_T200 = interpz(ds_hm['T'], ds_hm['z'], 200)
    hm_U200 = abs_velocity(hm_u200, hm_v200)

    # Read the Cabauw observations
    path = '/nobackup/users/stratum/Cabauw/'
    files = []
    for y in range(y0,y1+1):
        for m in range(m0,m1+1):
            files.append('{0:}/cesar_tower_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path,y,m))
    ds_cb = xr.open_mfdataset(files, drop_variables=['valid_dates'])

    # Calculate spectra
    freq = 1./600.

    cb_U_200 = Power_spectrum_lev(ds_cb['F'][:,0].values,  freq)
    hm_U_200 = Power_spectrum_lev(hm_U200,                 freq)
    cb_T_200 = Power_spectrum_lev(ds_cb['TA'][:,0].values, freq)
    hm_T_200 = Power_spectrum_lev(hm_T200,                 freq)


    pl.figure()
    pl.subplot(321)
    pl.loglog(cb_U_200.f[1:], cb_U_200.pxx[1:], linewidth=0.5, label='Cabauw')
    pl.loglog(hm_U_200.f[1:], hm_U_200.pxx[1:], linewidth=0.5, label='Harmonie')
    time_labels()
    pl.legend()

    pl.subplot(322)
    pl.loglog(cb_T_200.f[1:], cb_T_200.pxx[1:], linewidth=0.5, label='Cabauw')
    pl.loglog(hm_T_200.f[1:], hm_T_200.pxx[1:], linewidth=0.5, label='Harmonie')
    time_labels()
    pl.legend()

    pl.subplot(323)
    pl.loglog(cb_U_200.f[1:], cb_U_200.pxx_c[1:])
    pl.loglog(hm_U_200.f[1:], hm_U_200.pxx_c[1:])
    time_labels()

    pl.subplot(324)
    pl.loglog(cb_T_200.f[1:], cb_T_200.pxx_c[1:])
    pl.loglog(hm_T_200.f[1:], hm_T_200.pxx_c[1:])
    time_labels()

    pl.subplot(325)
    pl.semilogx(cb_U_200.f[1:], cb_U_200.pxx_c[1:] / cb_U_200.pxx_c[1])
    pl.semilogx(hm_U_200.f[1:], hm_U_200.pxx_c[1:] / hm_U_200.pxx_c[1])
    time_labels()

    pl.subplot(326)
    pl.semilogx(cb_T_200.f[1:], cb_T_200.pxx_c[1:] / cb_T_200.pxx_c[1])
    pl.semilogx(hm_T_200.f[1:], hm_T_200.pxx_c[1:] / hm_T_200.pxx_c[1])
    time_labels()
















