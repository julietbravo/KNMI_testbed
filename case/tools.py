import numpy as np
import datetime


def get_file_list(path, starttime, endtime):
    """
    Get list of required DDH NetCDF files to force
    LES from `starttime` to `endtime`
    """

    # For now limited to runs starting at a complete hour, to prevent having
    # to interpolate the inital conditions
    if starttime.minute != 0 or starttime.second != 0:
        raise RuntimeError('Can only create forcings starting at a complete hour!')

    # If experiment starts at start of cycle (t=0,3,6,..) we also need the previous cycle..
    if starttime.hour % 3 == 0:
        starttime -= datetime.timedelta(hours=3)

    # Number of cycles to convert
    n_cycles = int((endtime-starttime).total_seconds() / 3600. / 3.) + 1

    # Create list of cycles to convert
    files = []
    for i in range(n_cycles):
        date = starttime + i*datetime.timedelta(hours=3)
        in_file = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/LES_forcing_{1:04d}{2:02d}{3:02d}{4:02d}.nc'.\
            format(path, date.year, date.month, date.day, date.hour)
        files.append(in_file)

    return files


def get_start_end_indices(start, end, time):
    """
    Get indices in `time` that correspond to the requested `start` and `end` times
    """

    t0 = np.abs(np.datetime64(start) - time).argmin()
    t1 = np.abs(np.datetime64(end)   - time).argmin() + 1

    print('Using {} (index {}) to {} (index {})'.format(time[t0], t0, time[t1-1], t1-1))

    return t0, t1


def interpz(z_input, z_output, variable):
    """
    Interpolate (linear) `variable` from input grid (`z_input`) to output grid (`z_output`)
    """
    return np.interp(z_output, z_input, variable)


def interpz_time(z_input, z_output, variable):
    """
    Interpolate time varying `variable` from input grid (`z_input`) to output grid (`z_output`)
    """
    data = np.zeros((variable.shape[0], z_output.size))
    for t in range(variable.shape[0]):
        data[t,:] = interpz(z_input[t,:], z_output, variable[t,:])
    return data
