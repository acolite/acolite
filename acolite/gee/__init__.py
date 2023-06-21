from .agh import *
from .agh_old import *
from .agh_run import *
from .check_task import *
from .find_scenes import *

import ee
#ee.Authenticate() ## assume ee use is authenticated in current environment
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
