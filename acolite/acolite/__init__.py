from .acolite_gem import *

from . import settings

import os
path=os.path.dirname(os.path.abspath(__file__))
for i in range(2): path = os.path.split(path)[0]

## check if binary distribution
if '{}dist{}acolite'.format(os.path.sep,os.path.sep) in path:
    ## two more levels for this file
    for i in range(2): path = os.path.split(path)[0]
