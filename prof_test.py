# Run with the Profiler environment from the vscode debugger
# visualize with snakeviz: snakeviz /tmp/tmp.prof
import acolite as ac
import logging
from pathlib import Path
import cProfile, pstats
from datetime import datetime

now = datetime.now()
current_time = now.strftime("%H_%M_%S")

scene = "LC08_L1TP_114079_20220130_20220204_02_T1"
# Each scene will have all of its assets download into a sub-directory with it's id.
inputs_path = Path(f'/tmp/{scene}/inputs')
outputs_path = Path(f"/tmp/{scene}/outputs")
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