from .acolite_l1r import *
from .acolite_l2r import *
from .acolite_l2w import *

from .acolite_run import *
from .acolite_map import *
from .acolite_gui import *

from .acolite_luts import *

from .identify_bundle import *
from .inputfile_test import *
from .parameter_scaling import *

from . import settings
from . import logging

import os
path=os.path.dirname(os.path.abspath(__file__))
for i in range(2): path = os.path.split(path)[0]

## check if binary distribution
if '{}dist{}acolite'.format(os.path.sep,os.path.sep) in path:
    ## two more levels for this file
    for i in range(2): path = os.path.split(path)[0]
