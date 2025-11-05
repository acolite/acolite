from acolite import landsat
from acolite import sentinel2
from acolite import sentinel3
from acolite import planet
from acolite import pleiades
from acolite import worldview
from acolite import venus
from acolite import ikonos
from acolite import viirs
from acolite import seadas
from acolite import avhrr

from acolite import chris
from acolite import prisma
from acolite import hico
from acolite import hyperion
from acolite import desis
from acolite import enmap
from acolite import emit
from acolite import tanager
from acolite import hypso
from acolite import pace
from acolite import wyvern
from acolite import opencosmos
from acolite import earthcare

from acolite import seviri
from acolite import fci
from acolite import himawari
from acolite import goes
from acolite import goci


from acolite import gf
from acolite import amazonia
from acolite import formosat
from acolite import deimos

from acolite import ecostress
from acolite import haiyang
from acolite import huanjing
from acolite import sdgsat
from acolite import dimap
from acolite import s2resampling

from acolite import ac
from acolite import aerlut
from acolite import output
from acolite import shared
from acolite import dem
from acolite import ged
from acolite import masking

from acolite import tact
from acolite import acolite
from acolite import adjacency
from acolite import glint

from acolite import gem
from acolite import parameters

from acolite import api
from acolite import rtm
from acolite import custom
from acolite import zarr

## ignore numpy errors
import numpy as np
olderr = np.seterr(all='ignore')

## ignore osr exceptions
from osgeo import ogr
ogr.DontUseExceptions() #ogr.UseExceptions()

import os, sys, datetime, platform
## get platform identifiers
uname = platform.uname()
python = {'platform':sys.platform, 'version':sys.version}
system = {'sysname': uname.system, 'release': uname.release, 'machine': uname.machine, 'version': uname.version}

## find code_path and determine if binary
if getattr(sys, 'frozen', False): ## if pyinstaller binary
    binary = True #sys._MEIPASS
    code_path  = os.path.abspath(sys.argv[0]) ## use launch path
else:  ## if Python code
    binary = False
    code_path = os.path.dirname(os.path.abspath(__file__))
path = os.path.dirname(code_path)

## find config file
search_config = True
while search_config:
    cfile = os.path.join(path, 'config','config.txt')
    if os.path.exists(cfile):
        search_config = False
    else:
        path = os.path.split(path)[0]
    if len(path) <= 1:
        print('Could not find ACOLITE config file.')
        sys.exit()

## read config
config = shared.import_config(cfile)
config['code_path'] = code_path
config['path'] = path

## update version info
if 'version' in config:
    version = '{}'.format(config['version'])
else:
    version = 'GitHub Clone'

    gitdir = '{}/.git'.format(path)
    gd = {}
    if os.path.isdir(gitdir):
        gitfiles = os.listdir(gitdir)

        for f in ['ORIG_HEAD', 'FETCH_HEAD', 'HEAD']:
            gfile = '{}/{}'.format(gitdir, f)
            if not os.path.exists(gfile): continue
            st = os.stat(gfile)
            dt = datetime.datetime.fromtimestamp(st.st_mtime)
            gd[f] = dt.isoformat()[0:19]

        version_long = ''
        if 'HEAD' in gd:
            version_long+='clone {}'.format(gd['HEAD'])
            version = 'Generic GitHub Clone c{}'.format(gd['HEAD'])
        if 'FETCH_HEAD' in gd:
            version_long+=' pull {}'.format(gd['FETCH_HEAD'])
            version = 'Generic GitHub Clone p{}'.format(gd['FETCH_HEAD'])

## replace $ACDIR in config by ac.path
for t in config:
    if '$ACDIR' == config[t][0:6]:
        # os.path.join did not give the intended result on Windows
        config[t] = path + '/' + config[t].replace('$ACDIR', '')
        config[t] = config[t].replace('/', os.sep)
        ## make acolite dirs if they do not exist
        if not (os.path.exists(config[t])):
            os.makedirs(config[t])
## end replace $ACDIR

## add credentials
credentials = shared.import_config(config['credentials_file'])
for cr in credentials:
    if cr not in config: config[cr] = credentials[cr]
del credentials
## end add credentials

## run through config data
for t in config:
    ## split lists
    if ',' in config[t]:
        config[t] = config[t].split(',')
        continue
    ## set absolute paths
    if (os.path.exists(config[t])):
        config[t] = os.path.abspath(config[t])

## set attributes to skip
config['skip_attributes'] = acolite.settings.read_list(config['data_dir']+'/ACOLITE/attributes_skip.txt')

## set sensor lists
config['hyper_sensors'] = acolite.settings.read_list(config['data_dir']+'/ACOLITE/sensors_hyper.txt')
config['thermal_sensors'] = acolite.settings.read_list(config['data_dir']+'/ACOLITE/sensors_thermal.txt')
config['reverse_lut_sensors'] = acolite.settings.read_list(config['data_dir']+'/ACOLITE/sensors_reverse_lut.txt')

if 'verbosity' not in config: config['verbosity'] = 5

## read parameter scaling and settings
param = {'scaling':acolite.parameter_scaling(), 'discretisation': acolite.parameter_discretisation()}
import json
with open(config['parameter_cf_attributes'], 'r', encoding='utf-8') as f:
    param['attributes'] = json.load(f)

settings = {}
## read default processing settings
settings['defaults'] = acolite.settings.parse(None, settings=acolite.settings.load(None), merge=False)
## copy defaults to run, run will be updated with user settings and sensor defaults
settings['run'] = {k:settings['defaults'][k] for k in settings['defaults']}
## empty user settings
settings['user'] = {}
