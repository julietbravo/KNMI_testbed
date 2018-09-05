import pandas as pd
import glob
import datetime

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

        if start is not None and end is not None:
            df = df.loc[start:end]

    return df

if __name__ == '__main__':
    path = '/nobackup/users/stratum/FINO1_obs/'
    df   = read_FINO1(path)
