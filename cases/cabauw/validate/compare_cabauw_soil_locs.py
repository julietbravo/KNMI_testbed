import matplotlib.pyplot as pl
import pandas as pd
import xarray as xr
import numpy as np
import calendar
import datetime
import os

pl.close('all')

def read_Cabauw(y0, y1):

    # Read Cabauw observations
    pwd = '/nobackup/users/stratum/Cabauw'
    # List of files to read:
    files = []
    for year in range(y0, y1+1):
        for month in range(1,13):
            f = '{0:}/cesar_soil_water_lb1_t10_v1.1_{1:04d}{2:02d}.nc'.format(pwd, year, month)
            if os.path.exists(f):
                files.append(f)
            else:
                print('Can not find {}'.format(f))

    # Read with xarray:
    exclude = ['valid_dates', 'time_bnds', 'iso_dataset', 'product', 'station_details']
    cb = xr.open_mfdataset(files, concat_dim='time', drop_variables=exclude)

    # Convert to Pandas dataframe:
    cb = cb.to_dataframe()
    cb.index = cb.index.round('min')

    return cb

if __name__ == '__main__':
    if 'cb' not in locals():
        cb = read_Cabauw(2003,2017)

    cb_d = cb.groupby(cb.index.dayofyear).mean()

    if True:
        pl.figure()
        pl.subplot(131)
        for y in range(2003,2018):
            cb_y = cb.loc[cb.index.year == y]
            pl.plot(cb_y.index.dayofyear, cb_y['TH03'], '-', alpha=0.5)
        pl.plot(cb_d.index, cb_d['TH03'], 'k-', linewidth=2)
        pl.ylim(0.10, 0.9)

        pl.subplot(132)
        for y in range(2003,2018):
            cb_y = cb.loc[cb.index.year == y]
            pl.plot(cb_y.index.dayofyear, cb_y['TH08'], '-', alpha=0.5)
        pl.plot(cb_d.index, cb_d['TH08'], 'k-', linewidth=2)
        pl.ylim(0.10, 0.9)

        pl.subplot(133)
        for y in range(2003,2018):
            cb_y = cb.loc[cb.index.year == y]
            pl.plot(cb_y.index.dayofyear, cb_y['TH20'], '-', alpha=0.5)
        pl.plot(cb_d.index, cb_d['TH20'], 'k-', linewidth=2)
        pl.ylim(0.10, 0.9)


    if False:
        pl.figure()
        pl.subplot(121)
        pl.plot(cb.index, cb['TH03'], label='TH03 - 20030106-20170411')
        pl.plot(cb.index, cb['TH08'], label='TH08 - 20030106-20151212')
        pl.plot(cb.index, cb['TH20'], label='TH20 - 20030106-20110811')
        pl.legend()

        pl.subplot(122)
        for k,z in enumerate(z2):
            pl.plot(cb.index, cb['TH{0:02d}'.format(z)], label=str(z))
        pl.legend()






#    z1 = [3,8,20]
#    z2 = [5,19,33,40,56]


    #pl.figure()
    #pl.subplot(121)
    #pl.plot(cb.index, cb['TH03'], label='TH03 - 20030106-20170411')
    #pl.plot(cb.index, cb['TH08'], label='TH08 - 20030106-20151212')
    #pl.plot(cb.index, cb['TH20'], label='TH20 - 20030106-20170411')
    #pl.legend()
    #pl.xlim(datetime.datetime(2016,1,1), datetime.datetime(2018,1,1))

