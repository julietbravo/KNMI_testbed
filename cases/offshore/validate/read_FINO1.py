import matplotlib.pyplot as pl
import pandas as pd
import glob
import datetime
import numpy as np
from scipy import interpolate

class Correction_factors():
    def __init__(self):

        # Read correction factors from:
        # Westerhellweg A , Neumann T and Riedel V (2012) FINO1 Mast Correction. DEWI Magazin No. 40, February 2012
        low  = np.loadtxt('correction_factors_FINO1_low.txt', delimiter=',')
        high = np.loadtxt('correction_factors_FINO1_high.txt', delimiter=',')

        self.dir_low = low[:,0]
        self.fac_low = low[:,1]

        self.dir_high = high[:,0]
        self.fac_high = high[:,1]

        # Interpolation functions
        self.interp_low  = interpolate.interp1d(self.dir_low, self.fac_low)
        self.interp_high = interpolate.interp1d(self.dir_high, self.fac_high)

def wind_components(speed, direction):
    u = -speed * np.sin(np.deg2rad(direction))
    v = -speed * np.cos(np.deg2rad(direction))
    return u, v

def read_FINO1(path, start=None, end=None):

    files = glob.glob('{}/*.dat'.format(path))

    def add_vars(var_dict, name, shortname, heights):
        for z in heights:
            var_dict['{}_{}m'.format(name, z)] = '{}_{}'.format(shortname, z)

    # Dictionary for mapping long (German!) variable names to short ones
    variables = {}
    add_vars(variables, 'Globalstrahlung',      'Q',     [33,90])
    add_vars(variables, 'Luftdruck',            'p',     [20,90])
    add_vars(variables, 'Luftfeuchte',          'RH',    [33,40,50,70,90])
    add_vars(variables, 'Lufttemperatur',       'T',     [30,40,50,70,100])
    add_vars(variables, 'Windgeschwindigkeit',  'U',     [33,40,50,60,70,80,90,100])
    add_vars(variables, 'Windrichtung',         'Udir',  [33,40,50,60,70,80,90])

    for var in variables.keys():
        header = ['date', 'time', variables[var]]

        # Read data with Pandas
        tmp = pd.read_table('{}/FINO1_{}_20160101_20171231.dat'.format(path, var),\
                delim_whitespace=True, usecols=[0,1,2], names=header, skiprows=6,\
                parse_dates={'datetime' : ['date','time']}, na_values=[-999])
        tmp.set_index('datetime', inplace=True)

        # Create or merge data frames
        if 'df' not in locals():
            df = tmp
        else:
            df = df.merge(tmp, how='outer', left_index=True, right_index=True)

    # CHECK WITH INE!!!
    for z in [33,40,50,60,70,80,90]:
        var = df['Udir_{}'.format(z)]
        df['Udir_{}'.format(z)][var>360] -= 360

    # Calculate correction factors, and correct wind speeds
    corr = Correction_factors()
    for z in [40, 50, 60, 70]:
        df['fac_{}'.format(z)] = corr.interp_low(df['Udir_{}'.format(z)])
        df['Uc_{}' .format(z)] = df['U_{}' .format(z)] * df['fac_{}'.format(z)]
    for z in [80, 90]:
        df['fac_{}'.format(z)] = corr.interp_high(df['Udir_{}'.format(z)])
        df['Uc_{}' .format(z)] = df['U_{}' .format(z)] * df['fac_{}'.format(z)]

    # No correction functions available
    df['Uc_33' ] = df['U_33' ]
    df['Uc_100'] = df['U_100']

    # Corrected and uncorrected wind components
    for z in [33,40,50,60,70,80,90]:
        df['u_{}' .format(z)], df['v_{}' .format(z)] = wind_components(df['U_{}' .format(z)], df['Udir_{}'.format(z)])
        df['uc_{}'.format(z)], df['vc_{}'.format(z)] = wind_components(df['Uc_{}'.format(z)], df['Udir_{}'.format(z)])

    if start is not None and end is not None:
        df = df.loc[start:end]

    return df

if __name__ == '__main__':
    pl.close('all')

    #path = '/nobackup/users/stratum/FINO1_obs/'
    path = '/Users/bart/meteo/data/offshore_wind/FINO1/'

    if 'df' not in locals():
        df   = read_FINO1(path)

    # Mean corrected wind speed
    if (True):
        bins = np.arange(0, 360.01, 2)
        binc = 0.5 * (bins[1:] + bins[:-1])

        mean  = np.zeros_like(binc)
        meanc = np.zeros_like(binc)

        for i in range(binc.size):
            sel = ((df['Udir_90'] >= bins[i]) & (df['Udir_90'] < bins[i+1]))
            mean[i]  = df['U_90' ][sel].mean()
            meanc[i] = df['Uc_90'][sel].mean()

        pl.figure()
        ax=pl.subplot(polar=True)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        pl.plot(np.deg2rad(binc), mean,  label='uncorrected')
        pl.plot(np.deg2rad(binc), meanc, label='corrected')
        pl.legend()

    # Example correction factors and interpolation method
    if (False):
        correction = Correction_factors()

        pl.figure()
        pl.plot(correction.dir_low, correction.fac_low)

        dir = np.random.randint(low=0, high=360, size=20)
        fac = correction.interp_low(dir)

        pl.scatter(dir, fac)

