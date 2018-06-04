#!/usr/bin/env python

import glob
import os

def looks_like_python(name):
    if name.split('.')[-1] == 'py':
        return True
    else:
        return False

keep = ['namoptions.001', 'dales4']

files = glob.glob('*')

for f in files:
    if (f in keep or looks_like_python(f)):
        pass
    else:
        try:
            os.remove(f)
        except:
            print('Can not remove {}'.format(f))
