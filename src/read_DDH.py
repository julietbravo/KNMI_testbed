"""
  Quick 'n dirty Python "interface" to DDH (LFA) files.
  Should work for other LFA files, but `DDH_LFA`
  class contains some DDH specific routines.

  Requires LFA compiled, and its executables (lfac, lfalaf)
  available in the path. LFA code is available in external/LFA,
  and should compile out-of-the-box in KNMI desktops with the
  install script

  TO-DO: write real interface to LFA Fortran library!

  Bart van Stratum (KNMI)
"""


import subprocess
import numpy as np
import sys
import datetime
import os

# Extend PATH to include LFA binaries
path = os.path.abspath('{}/../external/LFA/'.format(os.path.dirname(os.path.abspath(__file__))))
os.environ["PATH"] += os.pathsep + path

# Some "private" functions
def _cl_call(call):
    sp = subprocess.Popen(call, shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
    return sp.stdout.read().decode("utf-8")


def _type_conversion(name):
    types = {'R4': np.float, 'I4': np.int, 'C': np.str}
    if name in types.keys():
        return types[name]
    else:
        print('Unknown type \"{}\", using np.float as default...'.format(name))
        return np.float


def _datetime_offset(date, offset, type):
    if type == 1:
        return date + offset*datetime.timedelta(hours=1)
    elif type == 2:
        return date + offset*datetime.timedelta(days=1)


class DDH_LFA:
    """
    Read individual DDH LFA files.
    Calls certain selected LFA command line routines,
    catches stdout output, and converts it to Python
    dictionaries or Numpy arrays.

    Arguments:
        file_path : path to DDH/LFA file
    """

    def __init__(self, file_path, levels=65):
        # Check if file actually exists
        if not os.path.exists(file_path):
            sys.exit('File {} does not exists...'.format(file_path))

        # Save path to file and file name
        self.file_path = file_path
        self.file_name = file_path.split('/')[-1]
        self.levels    = levels

        # Read the attributes
        self.read_attributes()


    def read_attributes(self):
        # Read `lfalaf` file dump
        data  = _cl_call('lfalaf {}'.format(self.file_path))
        lines = data.split('\n')[1:]

        # Loop over and process all records
        self.attributes = {}
        for line in lines:
            if (len(line) > 0):
                try:
                    bits   = line.replace('|','').split()
                    type_o = bits[1]                       # Original type
                    type   = _type_conversion(type_o)      # Numpy data type
                    length = int(bits[3])                  # Size of data array
                    name   = bits[4]                       # Variable name in LFA
                except: # not pretty to catch all, but LFA does weird stuff now and then..
                    print('Error reading LFA:')
                    print(data)
                    sys.exit()

                # Save variable as dictionary in dictionary....
                self.attributes[name] = {'type': type, 'length': length, 'type_o': type_o}

                # For `DOCFICHIER` and `DATE`, decode message
                if (name == 'DATE'):
                    date  = self.read_variable(name, reshape=False)
                    t_min = int(self.file_name.split('+')[-1])
                    start_date = datetime.datetime(year=date[0], month=date[1], day=date[2], hour=date[3])
                    valid_date = datetime.datetime(year=date[0], month=date[1], day=date[2], hour=date[3]) + datetime.timedelta(minutes=t_min)

                    self.attributes['datetime'] = {'start_date': start_date, 'forecast_date': valid_date}

                elif (name == 'DOCFICHIER'):
                    data  = self.read_variable(name, reshape=False)
                    self.attributes['doc'] = {'step': data[4], 'nlev': data[5]}


    def read_variable(self, name, reshape=True):
        if (name not in self.attributes.keys()):
            print('Variable \"{}\" does not exist in \"{}\".'.format(name, self.file_name))
            return None

        # Call `lfac` to get the data from one single variable:
        data = _cl_call('lfac {} {}'.format(self.file_path, name)).split('\n')

        # Remove empty elements
        cleaned_data = list(filter(None,data))

        # Return data as Numpy array of correct data type. `filter(None,data)` removes empty lines
        ncol = int(self.attributes[name]['length'] / self.levels)

        if len(cleaned_data) == 1:
            return np.array(cleaned_data, dtype=self.attributes[name]['type'])[0]
        elif (ncol%1 == 0 and reshape):
            return np.array(cleaned_data, dtype=self.attributes[name]['type']).reshape(int(ncol),-1).squeeze()
        else:
            return np.array(cleaned_data, dtype=self.attributes[name]['type'])


if __name__ == '__main__':
    # Only executed if script is called directly, for testing..

    data = '/Users/bart/meteo/data/Harmonie_DDH/20100228_00_ws/06/DHFDLHARM+0016'
    lfa  = DDH_LFA(data)
    v    = lfa.read_variable('VUU0')
