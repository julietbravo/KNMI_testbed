import xarray as xr

# --- Settings ---
f_in  = 'crossxzspan.004.nc'
t0    = 900
t1    = 1800
# ----------------

ds_in  = xr.open_dataset(f_in)
ds_out = ds_in.sel(time=slice(t0, t1)).mean(dim='time')
f_out  = f_in.split('.')[0] + '.mean.' + f_in.split('.')[1] + '.nc'
ds_out.to_netcdf(f_out)
