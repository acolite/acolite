# Run with the Profiler environment from the vscode debugger
# visualize with snakeviz: snakeviz /tmp/tmp.prof
import concurrent.futures
import cProfile
import logging
import pstats
import tempfile
from datetime import datetime
from pathlib import Path
from pprint import pprint
from time import time as timer
from urllib.parse import urlparse

import boto3
import geopandas as gpd
import pandas
import pystac_client
import shapely.geometry
from botocore.exceptions import ClientError
from threadpoolctl import threadpool_info
from tqdm.auto import tqdm

import acolite as ac

now = datetime.now()
current_time = now.strftime("%H_%M_%S")

scene = "LC08_L1TP_114079_20220130_20220204_02_T1"
# Each scene will have all of its assets download into a sub-directory with it's id.
inputs_path = Path(f'/tmp/sample_data/{scene}/inputs')
outputs_path = Path(f"/tmp/sample_data/{scene}/outputs")
inputs_path.mkdir(parents=True, exist_ok=True)

settings = {
    "inputfile": str(inputs_path),
    "output": str(outputs_path),
    # polygon=
    # limit=-29.9,152.5,-29.2,154.0
    "l2w_parameters": "Rrs_*",
    "rgb_rhot": False,
    "rgb_rhos": False,
    # "map_l2w": False,
    # "merge_zones": False,
    # "dsf_residual_glint_correction": True,
    "l2w_export_geotiff": True,
    "l1r_delete_netcdf" : True,
    "l2r_delete_netcdf": True,
    "l2t_delete_netcdf": True,
    "l2w_delete_netcdf": True,
    # "dsf_interface_reflectance": False, # False is the default
    "ancillary_data": False, # If you set this to True you must supply a username and password for EARTHDATA
    # "EARTHDATA_u": "",
    # "EARTHDATA_p": "",
    # "verbosity": 5
}

profiler = cProfile.Profile()
profiler.enable()

ac.acolite.acolite_run(settings=settings)

profiler.disable()
stats = pstats.Stats(profiler)
stats.dump_stats(f'./{current_time}_stats_file.dat')
stats.sort_stats('cumtime')
stats.strip_dirs()
stats.print_stats(15)