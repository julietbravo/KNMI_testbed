import matplotlib.pyplot as pl
import pandas as pd
import glob
import datetime
import numpy as np
from scipy import interpolate

class Correction_factors():
    def __init__(self):

        # Correction factors, read (manually....) from:
        # Westerhellweg A , Neumann T and Riedel V (2012) FINO1 Mast Correction. DEWI Magazin No. 40, February 2012
        x_40=np.array([  0.00,  18.04,  36.45,  55.61,  73.11,  89.18, 109.28, 129.39, 147.48, 163.87, 180.55, 197.82, 215.10, 236.85, 257.62, 273.61, 284.81, 293.72, 298.68, 301.99, 305.02, 307.89, 310.63, 312.06, 314.17, 316.27, 318.21, 319.33, 321.08, 323.74, 324.97, 330.05, 334.30, 336.25, 337.53, 339.48, 343.07, 347.23, 353.56, 360.00, ])
        y_40=np.array([0.9402, 0.9433, 0.9480, 0.9664, 0.9924, 1.0061, 1.0062, 1.0094, 1.0186, 1.0278, 1.0263, 1.0036, 0.9853, 0.9763, 0.9748, 0.9719, 0.9628, 0.9735, 1.0147, 1.0650, 1.1305, 1.2312, 1.3059, 1.3912, 1.4766, 1.5605, 1.6108, 1.6306, 1.6413, 1.6276, 1.6001, 1.5209, 1.4127, 1.3212, 1.2343, 1.1414, 1.0377, 0.9828, 0.9524, 0.9402, ])

        x_50=np.array([  0.00,  18.04,  36.45,  55.61,  73.11,  89.18, 109.28, 129.39, 147.48, 163.87, 180.55, 197.82, 215.10, 236.85, 257.62, 273.61, 284.81, 293.72, 298.68, 301.99, 305.02, 307.89, 310.63, 312.06, 314.17, 316.27, 318.21, 319.33, 321.08, 323.74, 324.97, 330.05, 334.30, 336.25, 337.53, 339.48, 343.07, 347.23, 353.56, 360.00, ])
        y_50=np.array([0.9402, 0.9433, 0.9480, 0.9664, 0.9924, 1.0061, 1.0062, 1.0094, 1.0186, 1.0278, 1.0263, 1.0036, 0.9853, 0.9763, 0.9748, 0.9719, 0.9628, 0.9735, 1.0147, 1.0650, 1.1305, 1.2312, 1.3059, 1.3912, 1.4766, 1.5605, 1.6108, 1.6306, 1.6413, 1.6276, 1.6001, 1.5209, 1.4127, 1.3212, 1.2343, 1.1414, 1.0377, 0.9828, 0.9524, 0.9402, ])

        x_60=np.array([  0.00,  18.04,  36.45,  55.61,  73.11,  89.18, 109.28, 129.39, 147.48, 163.87, 180.55, 197.82, 215.10, 236.85, 257.62, 273.61, 284.81, 293.72, 298.68, 301.99, 305.02, 307.89, 310.63, 312.06, 314.17, 316.27, 318.21, 319.33, 321.15, 324.89, 327.81, 329.94, 332.52, 334.30, 336.25, 337.53, 339.48, 343.07, 347.23, 353.56, 360.00, ])
        y_60=np.array([0.9402, 0.9433, 0.9480, 0.9664, 0.9924, 1.0061, 1.0062, 1.0094, 1.0186, 1.0278, 1.0263, 1.0036, 0.9853, 0.9763, 0.9748, 0.9719, 0.9628, 0.9735, 1.0147, 1.0650, 1.1305, 1.2312, 1.3059, 1.3912, 1.4766, 1.5605, 1.6108, 1.6306, 1.6550, 1.6550, 1.6230, 1.5712, 1.4691, 1.4127, 1.3212, 1.2343, 1.1414, 1.0377, 0.9828, 0.9524, 0.9402, ])

        x_70=np.array([  0.00,  18.04,  36.45,  55.61,  73.11,  89.18, 109.28, 129.39, 147.48, 163.87, 180.55, 197.82, 215.10, 236.85, 257.62, 273.61, 284.81, 293.72, 298.68, 301.99, 305.02, 307.89, 310.63, 312.06, 314.17, 316.27, 318.81, 323.21, 327.77, 330.59, 334.86, 336.49, 338.12, 340.43, 343.08, 347.03, 349.90, 355.94, 360.00, ])
        y_70=np.array([0.9402, 0.9433, 0.9480, 0.9664, 0.9924, 1.0061, 1.0062, 1.0094, 1.0186, 1.0278, 1.0263, 1.0036, 0.9853, 0.9763, 0.9748, 0.9719, 0.9628, 0.9735, 1.0147, 1.0650, 1.1305, 1.2312, 1.3059, 1.3912, 1.4766, 1.5605, 1.5925, 1.5864, 1.6154, 1.6352, 1.6017, 1.5163, 1.4279, 1.2709, 1.1109, 1.0118, 0.9707, 0.9524, 0.9433, ])

        x_80=np.array([  0.00,  16.68,  35.11,  55.98,  73.82,  92.31, 114.11, 140.33, 162.22, 186.31, 206.96, 226.65, 247.70, 267.06, 277.97, 287.24, 290.57, 293.73, 297.53, 301.08, 304.47, 307.73, 311.56, 315.92, 319.36, 322.55, 324.73, 326.91, 329.34, 331.86, 333.77, 335.58, 337.37, 340.71, 345.14, 352.18, 360.00, ])
        y_80=np.array([0.9402, 0.9433, 0.9510, 0.9725, 1.0000, 1.0199, 1.0200, 1.0201, 1.0385, 1.0203, 0.9929, 0.9793, 0.9657, 0.9535, 0.9566, 0.9719, 1.0283, 1.1183, 1.2022, 1.3028, 1.3714, 1.4126, 1.4293, 1.4156, 1.4218, 1.4492, 1.4782, 1.5072, 1.5163, 1.4736, 1.3730, 1.2526, 1.1276, 1.0438, 0.9722, 0.9493, 0.9402, ])

        x_90=np.array([  0.00,  18.40,  37.54,  52.96,  69.43,  84.51,  99.88, 118.25, 135.62, 153.39, 169.08, 187.37, 207.66, 224.94, 243.63, 258.53, 271.46, 280.38, 285.96, 288.25, 292.91, 295.04, 298.88, 301.55, 306.10, 313.50, 317.01, 319.60, 322.47, 324.57, 328.11, 330.38, 331.30, 333.78, 335.62, 339.06, 342.69, 351.50, 360.00, ])
        y_90=np.array([0.9494, 0.9479, 0.9617, 0.9816, 1.0061, 1.0275, 1.0352, 1.0307, 1.0323, 1.0445, 1.0476, 1.0264, 0.9960, 0.9778, 0.9687, 0.9520, 0.9490, 0.9627, 0.9887, 1.0405, 1.1610, 1.2509, 1.3424, 1.4019, 1.4278, 1.4080, 1.4294, 1.4736, 1.5041, 1.5163, 1.4721, 1.3776, 1.2846, 1.1611, 1.0468, 0.9828, 0.9584, 0.9493, 0.9463, ])

        self.interp_40 = interpolate.interp1d(x_40, y_40, kind='slinear')
        self.interp_50 = interpolate.interp1d(x_50, y_50, kind='slinear')
        self.interp_60 = interpolate.interp1d(x_60, y_60, kind='slinear')
        self.interp_70 = interpolate.interp1d(x_70, y_70, kind='slinear')
        self.interp_80 = interpolate.interp1d(x_80, y_80, kind='slinear')
        self.interp_90 = interpolate.interp1d(x_90, y_90, kind='slinear')


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
    for z in range(40,91,10):
        interp_func = getattr(corr, 'interp_{}'.format(z))

        df['fac_{}'.format(z)] = interp_func( df['Udir_{}'.format(z)] )             # Correction factor
        df['Uc_{}' .format(z)] = df['U_{}' .format(z)] * df['fac_{}'.format(z)]     # Corrected wind speed

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

    path = '/nobackup/users/stratum/FINO1_obs/'
    #path = '/Users/bart/meteo/data/offshore_wind/FINO1/'

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

