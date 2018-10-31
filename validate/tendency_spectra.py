import matplotlib.pyplot as pl
import xarray as xr
import pandas as pd
import numpy as np
from scipy.signal import periodogram

pl.ion()
pl.close('all')

# Path to the LES forcings
#path  = '/scratch/ms/nl/nkbs/DOWA/LES_forcing/'
#path  = '/nobackup/users/stratum/DOWA/LES_forcing/'
path   = '/Users/bart/meteo/data/Harmonie_LES_forcing/'

# Path to Cabauw observations
path_cb = '/Users/bart/meteo/observations/Cabauw'




def swin(doy, time, lat, lon):
    lon    = -lon
    sda    = 0.409 * np.cos(2. * np.pi * (doy - 173.) / 365.)
    sinlea = np.sin(2. * np.pi * lat / 360.) * np.sin(sda) - \
             np.cos(2. * np.pi * lat / 360.) * np.cos(sda) * \
             np.cos(2. * np.pi * (time*3600.) / 86400. - 2. * np.pi * lon / 360.)

    sinlea = np.maximum(sinlea, 1e-9)
    Tr     = (0.6 + 0.2 * sinlea)
    return 1368. * Tr * sinlea


def abs_velocity(u, v):
    return (u**2 + v**2)**0.5


def time_labels():
    ax = pl.gca()
    yext = ax.get_ybound()

    hours  = [1,3,24,24*30,24*365]
    labels = ['1h', '3h', 'day', 'month', 'year']

    for hour,label in zip(hours, labels):
        pl.plot([hour, hour], yext, 'C3', linewidth=2, alpha=0.4)
        pl.text(hour, yext[1], label, rotation=90, va='top', ha='right')

    ax.set_ylim(yext)


def read_harmonie(name):
    print('Reading {}'.format(name))

    # Read the merged (monthly) LES forcings
    files = []
    for y in range(2016,2018):
        for m in range(1,13):
            files.append('{0:}/{1:}_{2:04d}{3:02d}.nc'.format(path, name, y, m))

    return xr.open_mfdataset(files)


def read_cabauw(start_year, end_year):
    print('Reading Cabauw obs')

    files = []
    for y in range(start_year, end_year+1):
        print(y)
        for m in range(1,13):
            files.append('{0:}/cesar_tower_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path_cb, y, m))
            files.append('{0:}/cesar_surface_radiation_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path_cb, y, m))

    return xr.open_mfdataset(files, drop_variables=['valid_dates'], autoclose=True)


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
        ndims = len(array.shape)

        # Calculate the power spectrum (units = un**2, with `un` the units of `array`)
        self.k, self.pxx  = periodogram(array, fs=sampling_freq, scaling='spectrum', axis=0)

        # Vertical mean (TO-DO: height weighted)
        if ndims == 2:
            self.pxx_m = np.mean(self.pxx, axis=1)
        else:
            self.pxx_m = self.pxx

        # Cumulative
        self.pxx_c   = np.cumsum(self.pxx  [::-1], axis=0)[::-1]
        self.pxx_m_c = np.cumsum(self.pxx_m[::-1]        )[::-1]

        ## Frequency in hours:
        self.k[0] = 1e-16       # Mean == 0
        self.f = 1./self.k/3600.

        # Check-check-check....
        if ndims == 2:
            self.var_numpy = np.var(array[:,0])
            self.var_spec  = np.sum(self.pxx[1:,0])
        else:
            self.var_numpy = np.var(array)
            self.var_spec  = np.sum(self.pxx)
        print('Variance Numpy = {0:10.6e}, Periodogram = {1:10.6e}'.format(self.var_numpy, self.var_spec))
        #print('Variance Numpy = {0:10.6e}, Periodogram = {1:10.6e}, FFT = {2:10.6e}'.format(self.var_numpy, self.var_spec, var_fft))


class Power_spectra:
    def __init__(self, data, label, slice):
        self.label = '{}_{}_{}'.format(label, slice.start, slice.stop-1)
        start = slice.start
        stop  = slice.stop

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


if __name__ == '__main__':

    # ----------------------------
    # Spectra: Harmonie
    # ----------------------------
    if False:

        # Read the LES forcings
        if 'cb_30km' not in locals():
            cb_30km = read_harmonie('Cabauw_30km')
            lb_30km = read_harmonie('Loobos_30km')
            k1_30km = read_harmonie('K13_30km')
            f3_30km = read_harmonie('F3_30km')

        # Create power spectra
        if 'ps_cb' not in locals():
            slice = np.s_[0:1]       # 8=<200m, 21=<1000m

            ps_cb = Power_spectra(cb_30km, 'cb', slice)
            ps_lb = Power_spectra(lb_30km, 'lb', slice)
            ps_k1 = Power_spectra(k1_30km, 'k1', slice)
            ps_f3 = Power_spectra(f3_30km, 'f3', slice)


        # Selection of location & variable
        loc = ps_hm

        var1 = 'sU'; label1=r'$U$ (m s$^{-1}$)'; unit1=r'm$^2$ s$^{-2}$'
        var2 = 'sT'; label2=r'$T$ (K)';          unit2=r'K$^2$'

        #var1 = 'sdU'; label1=r'$\partial_t U$ (dynamics) (m s$^{-2}$)'; unit1=r'm$^2$ s$^{-4}$'
        #var2 = 'sdT'; label2=r'$\partial_t T$ (dynamics) (K s$^{-1}$)'; unit2=r'K$^2$ s$^{-2}$'

        #var1 = 'spU'; label1=r'$\partial_t U$ (physics) (m s$^{-2}$)'; unit1=r'm$^2$ s$^{-4}$'
        #var2 = 'spT'; label2=r'$\partial_t T$ (physics) (K s$^{-1}$)'; unit2=r'K$^2$ s$^{-2}$'

        #var1 = 'stU'; label1=r'$\partial_t U$ (total) (m s$^{-2}$)'; unit1=r'm$^2$ s$^{-4}$'
        #var2 = 'stT'; label2=r'$\partial_t T$ (total) (K s$^{-1}$)'; unit2=r'K$^2$ s$^{-2}$'

        pl.figure(figsize=(10,4))

        pl.subplot(121)
        pl.title(label1, loc='left')
        pl.loglog(loc.sU.f[1:], getattr(loc, var1).pxx_m[1:], linewidth=0.5, color='k')
        pl.ylabel(r'P ({})'.format(unit1))
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.subplot(122)
        pl.title(label2, loc='left')
        pl.loglog(loc.sU.f[1:], getattr(loc, var2).pxx_m[1:], linewidth=0.5, color='k')
        pl.ylabel(r'P ({})'.format(unit2))
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.tight_layout()
        pl.savefig('{}_P_{}_{}.pdf'.format(loc.label, var1, var2))

        pl.figure(figsize=(10,4))

        pl.subplot(121)
        pl.title(label1, loc='left')
        pl.loglog(loc.sU.f[1:], getattr(loc, var1).pxx_m_c[1:], linewidth=1.5, color='k')
        pl.ylabel(r'$\sum$ P ({})'.format(unit1))
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.subplot(122)
        pl.title(label2, loc='left')
        pl.loglog(loc.sU.f[1:], getattr(loc, var2).pxx_m_c[1:], linewidth=1.5, color='k')
        pl.ylabel(r'$\sum$ P ({})'.format(unit2))
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        time_labels()

        pl.tight_layout()
        pl.savefig('{}_cP_{}_{}.pdf'.format(loc.label, var1, var2))

        pl.figure(figsize=(10,4))

        pl.subplot(121)
        pl.title(label1, loc='left')
        pl.semilogx(loc.sU.f[1:], getattr(loc, var1).pxx_m_c[1:]/getattr(loc, var1).pxx_m_c[1], linewidth=1.5, color='k')
        pl.ylabel('$\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0,1)
        time_labels()

        pl.subplot(122)
        pl.title(label2, loc='left')
        pl.semilogx(loc.sU.f[1:], getattr(loc, var2).pxx_m_c[1:]/getattr(loc, var2).pxx_m_c[1], linewidth=1.5, color='k')
        pl.ylabel('$\sum$ P / total (-)')
        pl.xlabel('frequency (h)')
        pl.grid(axis='y', which='both')
        pl.ylim(0,1)
        time_labels()

        pl.tight_layout()
        pl.savefig('{}_cPn_{}_{}.pdf'.format(loc.label, var1, var2))




    # ----------------------------
    # Spectra: Harmonie vs Cabauw obs
    # ----------------------------
    if True:

        # Read Cabauw data
        if 'obs' not in locals():
            obs = read_cabauw(2001, 2017)

        # Read Harmonie data, and interpolate to obs heights
        if 'hm_data' not in locals():
            hm = read_harmonie('Cabauw_column')

            # Interpolate u, T, q to observation heights Cabauw, and create spectra
            hm_data = {}
            variables = ['T','u','v','U']
            for var in variables:
                print('Interpolating {}'.format(var))
                for z in obs.z.values:
                    if var == 'U':
                        hm_data['U_{0:03.0f}'.format(z)] = abs_velocity(hm_data['u_{0:03.0f}'.format(z)], hm_data['v_{0:03.0f}'.format(z)])
                    else:
                        hm_data['{0:}_{1:03.0f}'.format(var, z)] = interpz(hm[var], hm['z'], z)

        # Create spectra
        freq = 1./600.
        if 'hm_spec' not in locals():
            hm_spec = {}
            variables = ['T','U']
            for var in variables:
                for z in obs.z.values:
                    varname = '{0:}_{1:03.0f}'.format(var, z)
                    hm_spec[varname] = Power_spectrum(hm_data[varname], freq)

        if 'cb_spec' not in locals():
            cb_spec = {}
            variables = ['TA','F']
            for var in variables:
                for z in obs.z.values:
                    k = np.abs(obs.z.values - z).argmin()
                    varname = '{0:}_{1:03.0f}'.format(var, z)
                    cb_spec[varname] = Power_spectrum(obs[var][:,k].values, freq)

        if True:
            bins     = np.array([0,1,3,12,24,31*24,365*24])[::-1]
            index_hm = np.zeros_like(bins)
            index_cb = np.zeros_like(bins)

            labels = ['Y-M', 'M-D', '24-12h', '12-3h', '3-1h', '<1h']

            for i,freq in enumerate(bins):
                index_hm[i] = np.abs(hm_spec['U_010'].f - freq).argmin()
                index_cb[i] = np.abs(cb_spec['F_010'].f - freq).argmin()

            hm_var_T10 = np.zeros(bins.size-1)
            cb_var_T10 = np.zeros(bins.size-1)

            for i in range(bins.size-1):
                hm_var_T10[i] = np.sum(hm_spec['U_010'].pxx_m[index_hm[i]:index_hm[i+1]])
                cb_var_T10[i] = np.sum(cb_spec['F_010'].pxx_m[index_cb[i]:index_cb[i+1]])

            pl.figure()
            ax = pl.subplot(121)
            pl.bar(np.arange(bins.size-1)-0.2, hm_var_T10, width=0.3)
            pl.bar(np.arange(bins.size-1)+0.2, cb_var_T10, width=0.3)
            ax.set_xticklabels(labels)


        if False:
            pl.figure(figsize=(10,8))

            pl.subplot(221)
            pl.title('U 10 m', loc='left')
            pl.loglog(hm_spec['U_010'].f[1:], hm_spec['U_010'].pxx_m_c[1:], 'k-',  linewidth=1.5, label='Harmonie')
            pl.loglog(cb_spec['F_010'].f[1:], cb_spec['F_010'].pxx_m_c[1:], 'r--', linewidth=1.5, label='Observations')
            pl.legend()
            pl.ylabel(r'$\sum$ P (m$^2$ s$^{-2}$)')
            pl.xlabel('frequency (h)')
            pl.grid(True, axis='y', which='both')
            time_labels()

            pl.subplot(222)
            pl.title('T 10 m', loc='left')
            pl.loglog(hm_spec['T_010' ].f[1:], hm_spec['T_010' ].pxx_m_c[1:], 'k-',  linewidth=1.5)
            pl.loglog(cb_spec['TA_010'].f[1:], cb_spec['TA_010'].pxx_m_c[1:], 'r--', linewidth=1.5)
            pl.ylabel(r'$\sum$ P (K$^2$)')
            pl.xlabel('frequency (h)')
            pl.grid(True, axis='y', which='both')
            time_labels()

            pl.subplot(223)
            pl.title('U 200 m', loc='left')
            pl.loglog(hm_spec['U_200'].f[1:], hm_spec['U_200'].pxx_m_c[1:], 'k-',  linewidth=1.5)
            pl.loglog(cb_spec['F_200'].f[1:], cb_spec['F_200'].pxx_m_c[1:], 'r--', linewidth=1.5)
            pl.ylabel(r'$\sum$ P (m$^2$ s$^{-2}$)')
            pl.xlabel('frequency (h)')
            pl.grid(True, axis='y', which='both')
            time_labels()

            pl.subplot(224)
            pl.title('T 200 m', loc='left')
            pl.loglog(hm_spec['T_200' ].f[1:], hm_spec['T_200' ].pxx_m_c[1:], 'k-',  linewidth=1.5)
            pl.loglog(cb_spec['TA_200'].f[1:], cb_spec['TA_200'].pxx_m_c[1:], 'r--', linewidth=1.5)
            pl.ylabel(r'$\sum$ P (K$^2$)')
            pl.xlabel('frequency (h)')
            pl.grid(True, axis='y', which='both')
            time_labels()

            pl.tight_layout()

        if False:
            n = 1
            du_hm_10m = np.abs(hm_data['U_200'][n:]  - hm_data['U_200'][:-n])
            du_cb_10m = np.abs(obs['F'][n:,0].values - obs['F'][:-n,0].values)

            n = 2
            du_hm_20m = np.abs(hm_data['U_200'][n:]  - hm_data['U_200'][:-n])
            du_cb_20m = np.abs(obs['F'][n:,0].values - obs['F'][:-n,0].values)

            n = 3
            du_hm_30m = np.abs(hm_data['U_200'][n:]  - hm_data['U_200'][:-n])
            du_cb_30m = np.abs(obs['F'][n:,0].values - obs['F'][:-n,0].values)

            n = 6
            du_hm_60m = np.abs(hm_data['U_200'][n:]  - hm_data['U_200'][:-n])
            du_cb_60m = np.abs(obs['F'][n:,0].values - obs['F'][:-n,0].values)

            n = 12
            du_hm_120m = np.abs(hm_data['U_200'][n:]  - hm_data['U_200'][:-n])
            du_cb_120m = np.abs(obs['F'][n:,0].values - obs['F'][:-n,0].values)

            n = 18
            du_hm_180m = np.abs(hm_data['U_200'][n:]  - hm_data['U_200'][:-n])
            du_cb_180m = np.abs(obs['F'][n:,0].values - obs['F'][:-n,0].values)

            bins = np.arange(0,5.01,0.1)

            pl.figure(figsize=(10,6))
            pl.subplot(231)
            pl.title('10 min', loc='left', fontsize='small')
            pl.hist(du_hm_10m, bins, histtype='step', label='Harmonie')
            pl.hist(du_cb_10m, bins, histtype='step', label='Observations')
            pl.legend()
            pl.xlim(0,5)
            pl.xlabel(r'|$\Delta U$| (m s$^{-1}$)')
            pl.ylabel('N (-)')

            pl.subplot(232)
            pl.title('20 min', loc='left', fontsize='small')
            pl.hist(du_hm_20m, bins, histtype='step')
            pl.hist(du_cb_20m, bins, histtype='step')
            pl.xlim(0,5)
            pl.xlabel(r'|$\Delta U$| (m s$^{-1}$)')

            pl.subplot(233)
            pl.title('30 min', loc='left', fontsize='small')
            pl.hist(du_hm_30m, bins, histtype='step')
            pl.hist(du_cb_30m, bins, histtype='step')
            pl.xlim(0,5)
            pl.xlabel(r'|$\Delta U$| (m s$^{-1}$)')

            pl.subplot(234)
            pl.title('1 hour', loc='left', fontsize='small')
            pl.hist(du_hm_60m, bins, histtype='step')
            pl.hist(du_cb_60m, bins, histtype='step')
            pl.xlim(0,5)
            pl.xlabel(r'|$\Delta U$| (m s$^{-1}$)')
            pl.ylabel('N (-)')

            pl.subplot(235)
            pl.title('2 hours', loc='left', fontsize='small')
            pl.hist(du_hm_120m, bins, histtype='step')
            pl.hist(du_cb_120m, bins, histtype='step')
            pl.xlim(0,5)
            pl.xlabel(r'|$\Delta U$| (m s$^{-1}$)')

            pl.subplot(236)
            pl.title('3 hours', loc='left', fontsize='small')
            pl.hist(du_hm_180m, bins, histtype='step')
            pl.hist(du_cb_180m, bins, histtype='step')
            pl.xlim(0,5)
            pl.xlabel(r'|$\Delta U$| (m s$^{-1}$)')

            pl.tight_layout()

    # ----------------------
    # Fun with spectra... :)
    # ----------------------
    if False:

        def zero_except_list(fft, indices, zero_mean=False):
            fft = fft.copy()

            k0 = 0 if zero_mean else 1

            # Zero everything up to the first index
            fft.real[k0:indices[0]]   = 0
            fft.imag[k0:indices[0]]   = 0

            # Zero everything after the last index
            fft.real[indices[-1]+1:] = 0.
            fft.imag[indices[-1]+1:] = 0.

            # Zero everything in between
            for i in range(0,indices.size-1):
                fft.real[indices[i]+1:indices[i+1]] = 0.
                fft.imag[indices[i]+1:indices[i+1]] = 0.

            return fft


        def zero_except(fft, index):
            tmp = fft.copy()

            tmp.real[:index]   = 0.
            tmp.imag[:index]   = 0.
            tmp.real[index+1:] = 0.
            tmp.imag[index+1:] = 0.

            return tmp


        def spec_power(fft, n):
            fft = fft.copy() / n
            sp  = fft.real**2. + fft.imag**2.
            sp[1:n-1] *= 2.
            return sp


        class Decompose:
            def __init__(self, signal, wavenumbers):
                print('boe')

                self.n  = wavenumbers.size
                self.nt = signal.size

                # Signal from individual wavenumbers
                self.mode = np.zeros((self.n, self.nt), dtype=np.float)

                # 1. Fourier transform input signal
                fft = np.fft.rfft(signal)

                # Power spectra of original signal
                self.spf  = spec_power(fft, self.nt)
                self.k    = np.arange(int(self.nt/2+1))

                # Cumulative power spectra
                self.spfc = np.cumsum(self.spf[::-1])[::-1]

                # 2. Mean
                tmp = zero_except(fft, 0)
                self.mean = np.fft.irfft(tmp)[0]

                # Spectral power from individual wavenumbers
                self.sp = np.zeros((self.n, self.spf.size), dtype=np.float)

                # 3. Individual wavenumbers
                for i,wn in enumerate(wavenumbers):
                    tmp = zero_except(fft, wn)
                    self.mode[i,:] = np.fft.irfft(tmp)
                    self.sp[i,:] = spec_power(tmp, self.nt)

                # 4. Sum modes one-by-one
                self.modes = np.zeros_like(self.mode)
                self.modes[0,:] = self.mode[0,:] + self.mean
                for i in range(1,wavenumbers.size):
                    self.modes[i,:] = self.modes[i-1,:] + self.mode[i]


        # Read Cabauw data
        if 'obs' not in locals():
            obs   = read_cabauw(2001, 2017)
            n_obs = obs.dims['time']

        # Read Harmonie data, and interpolate to obs heights
        if 'cb_hm' not in locals():
            hm = read_harmonie('Cabauw_column')

        if 'T2' not in locals():
            # Wave numbers of peaks in spectra
            i24 = int(n_obs/24/6)    # f=24 h
            i12 = int(n_obs/12/6)    # f=12 h, etc.
            i8  = int(n_obs/8 /6)
            i6  = int(n_obs/6 /6)
            wns = np.array([i24, i12, i8, i6])

            # 2-meter temperature
            T2     = obs['TA'][:,-1].to_dataframe()
            T2mean = T2.groupby(T2.index.hour).mean()
            T2d    = Decompose(T2.values[:,1], wns)

            # 2-meter specific humidity
            Q2     = obs['Q'][:,-1].to_dataframe()
            Q2mean = Q2.groupby(Q2.index.hour).mean()
            Q2d    = Decompose(Q2.values[:,1], wns)

            # 10-meter wind
            F10     = obs['F'][:,-2].to_dataframe()
            F10mean = F10.groupby(F10.index.hour).mean()
            F10d    = Decompose(F10.values[:,1], wns)

            # Surface shortwave rad
            S     = obs['SWD'][:].to_dataframe()
            Smean = S.groupby(S.index.hour).mean()
            Sd    = Decompose(S.values[:,0], wns)

            # Time in hours
            time = np.arange(n_obs)/6.


        if True:
            ttot = float(obs.time[-1] - obs.time[0]) / 1e9 / 3600

            pl.figure(figsize=(10,7))

            pl.subplot(221)
            pl.title(r'a) $T_{2m}$', loc='left')
            pl.loglog(ttot/T2d.k, T2d.spf, linewidth=0.8, color='k', rasterized=True)
            time_labels()
            pl.xlabel(r'$f$ (h)')
            pl.ylabel(r'$P$ (K$^2$)')

            pl.subplot(222)
            pl.title(r'b) $q_{2m}$', loc='left')
            pl.loglog(ttot/Q2d.k, Q2d.spf, linewidth=0.8, color='k', rasterized=True)
            time_labels()
            pl.xlabel(r'$f$ (h)')
            pl.ylabel(r'$P$ (g$^2$ kg$^{-2}$)')

            pl.subplot(223)
            pl.title(r'c) $U_{10m}$', loc='left')
            pl.loglog(ttot/F10d.k, F10d.spf, linewidth=0.8, color='k', rasterized=True)
            time_labels()
            pl.xlabel(r'$f$ (h)')
            pl.ylabel(r'$P$ (m$^2$ s$^{-2}$)')

            pl.subplot(224)
            pl.title(r'd) $SW_\mathrm{in}$', loc='left')
            pl.loglog(ttot/Sd.k, Sd.spf, linewidth=0.8, color='k', rasterized=True)
            time_labels()
            pl.xlabel(r'$f$ (h)')
            pl.ylabel(r'$P$ (W$^2$ m$^{-4}$)')

            pl.tight_layout()
            pl.savefig('P_Cabauw_obs.pdf')


            pl.figure(figsize=(10,7))

            pl.subplot(221)
            pl.title(r'a) $T_{2m}$', loc='left')
            pl.semilogx(ttot/T2d.k, T2d.spfc/T2d.spfc[1], color='k', rasterized=True)
            time_labels()
            pl.xlabel(r'$f$ (h)')
            pl.ylabel(r'$\hat{P}$ (-)')

            pl.subplot(222)
            pl.title(r'b) $q_{2m}$', loc='left')
            pl.semilogx(ttot/Q2d.k, Q2d.spfc/Q2d.spfc[1], color='k', rasterized=True)
            time_labels()
            pl.xlabel(r'$f$ (h)')
            pl.ylabel(r'$\hat{P}$ (-)')

            pl.subplot(223)
            pl.title(r'c) $U_{10m}$', loc='left')
            pl.semilogx(ttot/F10d.k, F10d.spfc/F10d.spfc[1], color='k', rasterized=True)
            time_labels()
            pl.xlabel(r'$f$ (h)')
            pl.ylabel(r'$\hat{P}$ (-)')

            pl.subplot(224)
            pl.title(r'd) $SW_\mathrm{in}$', loc='left')
            pl.semilogx(ttot/Sd.k, Sd.spfc/Sd.spfc[1], color='k', rasterized=True)
            time_labels()
            pl.xlabel(r'$f$ (h)')
            pl.ylabel(r'$\hat{P}$ (-)')

            pl.tight_layout()
            pl.savefig('Pn_Cabauw_obs.pdf')


        if False:
            day = np.s_[0:24*6]

            pl.figure(figsize=(10,4))

            pl.subplot(121)
            pl.title(r'$\mu$ = {0:.2f} K'.format(T2d.mean), loc='left', fontsize='small')
            pl.plot(T2mean.index+0.5, T2mean.TA, 'ks', label='Mean daily cycle', linewidth=2, mfc='none')
            pl.plot(time[day], T2d.modes[0,day], label=r'$\mu$ + $f$=24 h')
            pl.plot(time[day], T2d.modes[1,day], label=r'$\mu$ + $f$=24+12 h')
            pl.plot(time[day], T2d.modes[2,day], label=r'$\mu$ + $f$=24+12+8 h')
            pl.xlim(0,24)
            pl.xticks([0,6,12,18,24])
            pl.legend(loc='upper left', fontsize='small')
            pl.xlabel('time UTC')
            pl.ylabel(r'$T_\mathrm{2m}$ (K)')

            pl.subplot(122)
            pl.plot(time[day], T2d.mode[0,day], label=r'$f$=24 h')
            pl.plot(time[day], T2d.mode[1,day], label=r'$f$=12 h')
            pl.plot(time[day], T2d.mode[2,day], label=r'$f$=8 h')
            pl.plot(time[day], T2d.mode[3,day], label=r'$f$=6 h')
            pl.xlim(0,24)
            pl.xticks([0,6,12,18,24])
            pl.legend(fontsize='small')
            pl.xlabel('time UTC')
            pl.ylabel(r'$\Delta T_\mathrm{2m}$ (K)')

            pl.tight_layout()
            pl.savefig('modes_T2m_Cabauw.pdf')


            pl.figure(figsize=(10,4))

            pl.subplot(121)
            pl.title(r'$\mu$ = {0:.2f} g kg$^{{-1}}$'.format(Q2d.mean), loc='left', fontsize='small')
            pl.plot(Q2mean.index+0.5, Q2mean.Q, 'ks', label='Mean daily cycle', linewidth=2, mfc='none')
            pl.plot(time[day], Q2d.modes[0,day], label=r'$\mu$ + $f$=24 h')
            pl.plot(time[day], Q2d.modes[1,day], label=r'$\mu$ + $f$=24+12 h')
            pl.plot(time[day], Q2d.modes[2,day], label=r'$\mu$ + $f$=24+12+8 h')
            pl.xlim(0,24)
            pl.xticks([0,6,12,18,24])
            pl.legend(loc='upper left', fontsize='small')
            pl.xlabel('time UTC')
            pl.ylabel(r'$q_\mathrm{2m}$ (g kg$^{-1}$)')

            pl.subplot(122)
            pl.plot(time[day], Q2d.mode[0,day], label=r'$f$=24 h')
            pl.plot(time[day], Q2d.mode[1,day], label=r'$f$=12 h')
            pl.plot(time[day], Q2d.mode[2,day], label=r'$f$=8 h')
            pl.plot(time[day], Q2d.mode[3,day], label=r'$f$=6 h')
            pl.xlim(0,24)
            pl.xticks([0,6,12,18,24])
            pl.legend(fontsize='small')
            pl.xlabel('time UTC')
            pl.ylabel(r'$\Delta q_\mathrm{2m}$ (g kg$^{-1}$)')

            pl.tight_layout()
            pl.savefig('modes_Q2m_Cabauw.pdf')


            pl.figure(figsize=(10,4))

            pl.subplot(121)
            pl.title(r'$\mu$ = {0:.2f} m s$^{{-1}}$'.format(F10d.mean), loc='left', fontsize='small')
            pl.plot(F10mean.index+0.5, F10mean.F, 'ks', label='Mean daily cycle', linewidth=2, mfc='none')
            pl.plot(time[day], F10d.modes[0,day], label=r'$\mu$ + $f$=24 h')
            pl.plot(time[day], F10d.modes[1,day], label=r'$\mu$ + $f$=24+12 h')
            pl.plot(time[day], F10d.modes[2,day], label=r'$\mu$ + $f$=24+12+8 h')
            pl.xlim(0,24)
            pl.xticks([0,6,12,18,24])
            pl.legend(loc='upper left', fontsize='small')
            pl.xlabel('time UTC')
            pl.ylabel(r'$U_\mathrm{10m}$ (m/s)')

            pl.subplot(122)
            pl.plot(time[day], F10d.mode[0,day], label=r'$f$=24 h')
            pl.plot(time[day], F10d.mode[1,day], label=r'$f$=12 h')
            pl.plot(time[day], F10d.mode[2,day], label=r'$f$=8 h')
            pl.plot(time[day], F10d.mode[3,day], label=r'$f$=6 h')
            pl.xlim(0,24)
            pl.xticks([0,6,12,18,24])
            pl.legend(fontsize='small')
            pl.xlabel('time UTC')
            pl.ylabel(r'$\Delta U_\mathrm{10m}$ (m/s)')

            pl.tight_layout()
            pl.savefig('modes_U10m_Cabauw.pdf')


            pl.figure(figsize=(10,4))

            time_sw = np.linspace(0,24,256)
            sw1 = swin(172, time_sw, 51.97, 4.9)
            sw2 = swin(356, time_sw, 51.97, 4.9)

            pl.subplot(121)
            pl.title(r'$\mu$ = {0:.2f} W m$^{{-2}}$'.format(Sd.mean), loc='left', fontsize='small')
            pl.plot(Smean.index+0.5, Smean.SWD, 'ks', label='Mean daily cycle', linewidth=2, mfc='none')
            pl.plot(time[day], Sd.modes[0,day], label=r'$\mu$ + $f$=24 h')
            pl.plot(time[day], Sd.modes[1,day], label=r'$\mu$ + $f$=24+12 h')
            pl.plot(time[day], Sd.modes[2,day], label=r'$\mu$ + $f$=24+12+8 h')
            pl.xlim(0,24)
            #pl.ylim(-80,400)
            pl.xticks([0,6,12,18,24])
            pl.legend(loc='upper left', fontsize='small')
            pl.xlabel('time UTC')
            pl.ylabel(r'$\mathrm{SW}_\mathrm{in}$ (W/m2)')

            pl.subplot(122)
            pl.plot(time[day], Sd.mode[0,day], label=r'$f$=24 h')
            pl.plot(time[day], Sd.mode[1,day], label=r'$f$=12 h')
            pl.plot(time[day], Sd.mode[2,day], label=r'$f$=8 h')
            pl.plot(time[day], Sd.mode[3,day], label=r'$f$=6 h')
            pl.xlim(0,24)
            pl.xticks([0,6,12,18,24])
            pl.legend(fontsize='small')
            pl.xlabel('time UTC')
            pl.ylabel(r'$\Delta \mathrm{SW}_\mathrm{in}$ (W/m2)')

            pl.tight_layout()
            pl.savefig('modes_SWin_Cabauw.pdf')





        #time_sw = np.linspace(0,24,2048)
        #daylen = np.zeros(365)
        #for d in range(1,366):
        #    sw = swin(d, time_sw, 51.97, 4.9)

        #    dsw = sw[1:] - sw[:-1]
        #    pos = np.where(np.abs(dsw) > 1)[0]
        #    t0  = pos[0]
        #    t1  = pos[-1]

        #    daylen[d-1] = time_sw[t1] - time_sw[t0]

        #print(daylen.mean())




        # Temperature
        # --------------
        #fft = np.fft.rfft(T10)

        ## Including mean
        #fft1 = zero_except(fft, np.array([i24]))
        #fft2 = zero_except(fft, np.array([i24, i12]))
        #fft3 = zero_except(fft, np.array([i24, i12, i8]))
        #fft4 = zero_except(fft, np.array([i24, i12, i8, i6]))

        #T101 = np.fft.irfft(fft1)
        #T102 = np.fft.irfft(fft2)
        #T103 = np.fft.irfft(fft3)
        #T104 = np.fft.irfft(fft4)

        ## Excluding mean
        #fft11 = zero_except(fft, np.array([i24]), True)
        #fft21 = zero_except(fft, np.array([i12]), True)
        #fft31 = zero_except(fft, np.array([i8]),  True)
        #fft41 = zero_except(fft, np.array([i6]),  True)

        #T1011 = np.fft.irfft(fft11)
        #T1021 = np.fft.irfft(fft21)
        #T1031 = np.fft.irfft(fft31)
        #T1041 = np.fft.irfft(fft41)

        ## Wind
        #fft = np.fft.rfft(F10)

        ## Including mean
        #fft1 = zero_except(fft, np.array([i24]))
        #fft2 = zero_except(fft, np.array([i24, i12]))
        #fft3 = zero_except(fft, np.array([i24, i12, i8]))
        #fft4 = zero_except(fft, np.array([i24, i12, i8, i6]))

        #F101 = np.fft.irfft(fft1)
        #F102 = np.fft.irfft(fft2)
        #F103 = np.fft.irfft(fft3)
        #F104 = np.fft.irfft(fft4)

        ## Exluding mean
        #fft11 = zero_except(fft, np.array([i24]), True)
        #fft21 = zero_except(fft, np.array([i12]), True)
        #fft31 = zero_except(fft, np.array([i8]),  True)
        #fft41 = zero_except(fft, np.array([i6]),  True)

        #F1011 = np.fft.irfft(fft11)
        #F1021 = np.fft.irfft(fft21)
        #F1031 = np.fft.irfft(fft31)
        #F1041 = np.fft.irfft(fft41)

        ## Radiation
        #fft = np.fft.rfft(Q10)

        #fft1 = zero_except(fft, np.array([i24]))
        #fft2 = zero_except(fft, np.array([i24, i12]))
        #fft3 = zero_except(fft, np.array([i24, i12, i8]))
        #fft4 = zero_except(fft, np.array([i24, i12, i8, i6]))

        #Q101 = np.fft.irfft(fft1)
        #Q102 = np.fft.irfft(fft2)
        #Q103 = np.fft.irfft(fft3)
        #Q104 = np.fft.irfft(fft4)

        ## Exluding mean
        #fft11 = zero_except(fft, np.array([i24]), True)
        #fft21 = zero_except(fft, np.array([i12]), True)
        #fft31 = zero_except(fft, np.array([i8]),  True)
        #fft41 = zero_except(fft, np.array([i6]),  True)

        #Q1011 = np.fft.irfft(fft11)
        #Q1021 = np.fft.irfft(fft21)
        #Q1031 = np.fft.irfft(fft31)
        #Q1041 = np.fft.irfft(fft41)




        """
        pl.figure(figsize=(10,8))
        pl.subplot(221)
        pl.plot(Tmean.index+0.5, Tmean.TA, 'ks', label='Mean signal', linewidth=2, mfc='none')
        pl.plot(time, T101, '-', label='mean + freq=24h')
        pl.plot(time, T102, '-', label='mean + freq=24+12h')
        pl.plot(time, T103, '-', label='mean + freq=24+12+8h')
        pl.xlim(0,24)
        pl.legend()
        pl.xlabel('time UTC')
        pl.ylabel('T 2m (K)')

        pl.subplot(222)
        pl.plot(time, T1011, '-', label='freq=24h')
        pl.plot(time, T1021, '-', label='freq=12h')
        pl.plot(time, T1031, '-', label='freq=8h')
        pl.plot(time, T1041, '-', label='freq=6h')
        pl.xlim(0,24)
        pl.legend()
        pl.xlabel('time UTC')
        pl.ylabel('dT 2m (K)')

        pl.subplot(223)
        pl.plot(Fmean.index+0.5, Fmean.F, 'ks', label='Mean signal', linewidth=2, mfc='none')
        pl.plot(time, F101, '-', label='mean + freq=24h')
        pl.plot(time, F102, '-', label='mean + freq=24+12h')
        pl.plot(time, F103, '-', label='mean + freq=24+12+8h')
        pl.xlim(0,24)
        pl.legend()
        pl.xlabel('time UTC')
        pl.ylabel('U 10m (m/s)')

        pl.subplot(224)
        pl.plot(time, F1011, '-', label='freq=24h')
        pl.plot(time, F1021, '-', label='freq=12h')
        pl.plot(time, F1031, '-', label='freq=8h')
        pl.plot(time, F1041, '-', label='freq=6h')
        pl.xlim(0,24)
        pl.legend()
        pl.xlabel('time UTC')
        pl.ylabel('dU 10m (m/s)')

        pl.tight_layout()




        pl.figure(figsize=(10,8))
        pl.subplot(121)
        pl.plot(Qmean.index+0.5, Qmean.SWD, 'ks', label='Mean signal', linewidth=2, mfc='none')
        pl.plot(time, Q101, '-', label='mean + freq=24h')
        pl.plot(time, Q102, '-', label='mean + freq=24+12h')
        pl.plot(time, Q103, '-', label='mean + freq=24+12+8h')
        pl.xlim(0,24)
        pl.legend()
        pl.xlabel('time UTC')
        pl.ylabel('SWD (W m-2)')

        pl.subplot(122)
        pl.plot(time, Q1011, '-', label='freq=24h')
        pl.plot(time, Q1021, '-', label='freq=12h')
        pl.plot(time, Q1031, '-', label='freq=8h')
        pl.plot(time, Q1041, '-', label='freq=6h')
        pl.xlim(0,24)
        pl.legend()
        pl.xlabel('time UTC')
        pl.ylabel('dSWD (W m-2)')

        """


    # ------------- Supporting figures -------------
    if False:

        #case = cb_30km
        #name = 'cb'

        case = f3_30km
        name = 'f3'

        levs = np.array([4,14,20])
        z    = np.array([100,500,1000])

        # Montly averaged dynamic temperature tendency
        dTm  = np.zeros((levs.size, 12))
        for k,lev in enumerate(levs):
            df = case['dtT_dyn'].sel(level=lev).to_dataframe()
            dTm[k,:] = np.squeeze( df.groupby(df.index.month).mean() )

        # Daily cycle dynamic temperature tendency
        dTd_mjja = np.zeros((levs.size, 24))
        dTd_ndjf = np.zeros((levs.size, 24))
        for k,lev in enumerate(levs):
            df   = case['dtT_dyn'].sel(level=lev).to_dataframe()

            df_s = df.loc[(df.index.month >= 5)&(df.index.month <=8)]
            df_w = df.loc[(df.index.month >= 11)|(df.index.month <=2)]

            dTd_mjja[k,:] = np.squeeze( df_s.groupby(df_s.index.hour).mean() )
            dTd_ndjf[k,:] = np.squeeze( df_w.groupby(df_w.index.hour).mean() )


        pl.figure(figsize=(10,4))

        ax=pl.subplot(121)
        for k in range(levs.size):
            pl.plot(np.arange(1,13), dTm[k,:]*3600, label='z={} m'.format(z[k]))
        pl.legend()
        pl.ylabel(r'$\partial_t T$ (dynamics) (K h$^{-1}$)')
        pl.xlabel('month (-)')
        pl.xlim(1,12)
        ax.set_xticks(np.arange(1,13))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])

        pl.subplot(122)
        pl.plot(np.arange(0,24), dTd_mjja[0,:]*3600, label='z=100 m, MJJA')
        pl.plot(np.arange(0,24), dTd_ndjf[0,:]*3600, label='z=100 m, NDJF')
        pl.legend()
        pl.ylabel(r'$\partial_t T$ (dynamics) (K h$^{-1}$)')
        pl.xlabel('time UTC (h)')
        pl.xlim(0,23)

        pl.tight_layout()
        pl.savefig('{}_budget.pdf'.format(loc.label))


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
















