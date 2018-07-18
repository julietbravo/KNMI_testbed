import matplotlib.pyplot as pl
import pandas as pd
import numpy as np

from read_DDH_netcdf import *

def read_LEG_EPL(path):
    # Read and clean header
    f     = open(path, 'r', encoding="ISO-8859-1")
    info  = f.readline(). split(';')
    names = f.readline().split(';')
    f.close()
    
    for i in range(len(names)):
        names[i] = names[i].replace('LEG_', '')
        names[i] = names[i].replace('EPL_', '')
    
    return pd.read_csv(path, skiprows=3, header=1, sep=';', index_col=0, names=names, parse_dates=True, na_values=9999)

def mask_nan(array):
    return np.ma.masked_invalid(array)

def abs_speed(u,v):
    return (u**2 + v**2)**0.5


if __name__ == '__main__':

    if 'data_is_loaded' not in locals():
        data_is_loaded=True
        
        # Period to study
        start = datetime.datetime(year=2017, month=1, day=1,  hour=0)
        end   = datetime.datetime(year=2017, month=1, day=14, hour=21)

        # Harmonie data from DDH files
        # ----------------------------
        path = '/Users/bart/meteo/data/Harmonie_DDH/'
        ham  = read_DDH_netcdf(start, end, path)

        # Lichteiland Goeree observations
        # -------------------------------
        leg = read_LEG_EPL('/Users/bart/meteo/data/offshore_wind/LEG/statistics_2017.csv')
        leg = leg.loc[start:end]

        # Lichteiland Goeree observations
        # -------------------------------
        epl = read_LEG_EPL('/Users/bart/meteo/data/offshore_wind/EPL/statistics_2017.csv')
        epl = epl.loc[start:end]



    iloc_leg = 1
    iloc_epl = 2
    k = 2

    pl.close('all')
    pl.figure()
    pl.subplot(121)
    pl.plot(leg.index, mask_nan(leg['H063_Ws']), label='LEG')
    pl.plot(ham.time,  abs_speed(ham['u'][:,iloc_leg,k], ham['v'][:,iloc_epl,k]), label='Harmonie')
    pl.legend()

    pl.subplot(122)
    pl.plot(epl.index, mask_nan(epl['H063_WsHor_avg']), label='EPL')
    pl.plot(ham.time,  abs_speed(ham['u'][:,iloc_epl,k], ham['v'][:,iloc_epl,k]), label='Harmonie')
    pl.legend()
