#!/usr/bin/env python

import glob
import os

def looks_like_python(name):
    if name.split('.')[-1] == 'py':
        return True
    else:
        return False

keep = ['namoptions.001','namoptions.002','dales4','run.PBS','rrtmg_sw.nc','rrtmg_lw.nc']

files = glob.glob('*')

for f in files:
    if not (f in keep or looks_like_python(f)):
        try:
            os.remove(f)
        except:
            print('Can not remove {}'.format(f))
