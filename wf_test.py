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


def fetch_s3_obj(entry):
    path, uri, requester_pays = entry
    uu = urlparse(uri)
    bucket = uu.netloc
    key = uu.path.lstrip("/")
    requestPayer = 'requester' if requester_pays else ''
    session = boto3.session.Session() # Create a new session for each thread for thread-safety
    s3_client = session.client('s3')

    s3_client.download_file(bucket,
                            key,
                            str(path),
                            ExtraArgs = {'RequestPayer': requestPayer})
    return path


def s3_upload_file(file_name, bucket, object_name=None):
    """Upload a file to an S3 bucket

    :param file_name: File to upload
    :param bucket: Bucket to upload to
    :param object_name: S3 object name. If not specified then file_name is used
    :return: True if file was uploaded, else False
    """

    # If S3 object_name was not specified, use file_name
    if object_name is None:
        object_name = file_name

    # boto3
    if 's3' not in locals():
        s3 = boto3.client('s3')
    # Upload the file
    try:
        s3.upload_file(file_name, bucket, object_name)  # Config=config
    except (ClientError, FileNotFoundError) as e:
        logging.error(e)
        return False
    return True

def convert_bounds(bbox, invert_y=False):
    """
    Helper method for changing bounding box representation to leaflet notation

    ``(lon1, lat1, lon2, lat2) -> ((lat1, lon1), (lat2, lon2))``
    """
    x1, y1, x2, y2 = bbox
    if invert_y:
        y1, y2 = y2, y1
    return ((y1, x1), (y2, x2))


def download_and_process(stac_item, userid, scratch_bucket, s3_prefix):
    logger = logging.getLogger("mylogger")
    logger.setLevel(logging.DEBUG)
    if len(logger.handlers) == 0:
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        logger.addHandler(ch)
    logger.debug("download_and_process")

    # Each scene will have all of its assets download into a sub-directory with it's id.
    with tempfile.TemporaryDirectory(prefix=stac_item['id']) as inputs_path:
        with tempfile.TemporaryDirectory(prefix=stac_item['id']) as outputs_path:
            try:
                # Filter assets to be that of data (no overviews of html pages design for humans)
                assets = stac_item["assets"]
                assets_df = pandas.DataFrame.from_dict(assets, orient="index")
                filtered_assets = assets_df[
                    (assets_df["type"] != "text/html")
                    & (assets_df["type"] != "image/jpeg")
                ]["alternate"].to_list()

                # Create a list of things to download and process on Dask Workers
                # s3_item is a list of tuples (filepath, url, requester_pays)
                s3_items = [
                    (
                        Path(
                            inputs_path, Path(urlparse(s3_obj["s3"]["href"]).path).name
                        ),
                        s3_obj["s3"]["href"],
                        s3_obj["s3"]["storage:requester_pays"],
                    )
                    for s3_obj in filtered_assets
                ]

                # logger.debug("download items")
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    for p in executor.map(fetch_s3_obj, s3_items):
                        # logger.debug(p)
                        pass

                logger.debug("run acolite")

                # run Acolite
                settings = {
                    "inputfile": inputs_path,
                    "output": outputs_path,
                    # polygon=
                    # limit=-29.9,152.5,-29.2,154.0
                    "l2w_parameters": "Rrs_*",
                    "l2w_data_in_memory": True,
                    "rgb_rhot": False,
                    "rgb_rhos": False,
                    "map_l2w": False,
                    "merge_zones": False,
                    "dsf_residual_glint_correction": True,
                    "l2w_export_geotiff": True,
                    "l1r_delete_netcdf": True,
                    "l2r_delete_netcdf": True,
                    "l2t_delete_netcdf": True,
                    "l2w_delete_netcdf": True,
                    "dsf_interface_reflectance": False,  # False is the default
                    "ancillary_data": False,  # If you set this to True you must supply a username and password for EARTHDATA
                    "EARTHDATA_u": "",
                    "EARTHDATA_p": "",
                    "verbosity": 5,
                }
                ac.acolite.acolite_run(settings=settings)

                # Store outputs in user scratch
                # root = Path(settings["output"])
                # for file in root.iterdir():
                #     logging.debug(str(file.relative_to(root)))
                #     res = s3_upload_file(
                #         str(file),
                #         scratch_bucket,
                #         f'{userid}/{s3_outputs_prefix}/{stac_item["id"]}/{str(file.relative_to(root))}',
                #     )
                #     if not res:
                #         raise Exception(f'FAILED: {stac_item["id"]}')
            except (FileExistsError) as e:
                return f"{e}"
            except Exception as e:
                return f"{e}"

    return f's3://{scratch_bucket}/{userid}/{s3_outputs_prefix}/{stac_item["id"]}'









now = datetime.now()
current_time = now.strftime("%H_%M_%S")


x, y = (146.061975, -42.815608)  # Center point of a query
km2deg = 1.0 / 111
r = 25 * km2deg
bbox = (x - r, y - r, x + r, y + r)


catalog = pystac_client.Client.open("https://landsatlook.usgs.gov/stac-server")


query = catalog.search(
    collections=['landsat-c2l1'], datetime="2022-01-01/2022-01-31", max_items=None, limit=1, bbox=bbox
)

for ic in query.pages():
    # PySTAC ItemCollection
    # Dictionary (GeoJSON FeatureCollection)
    items_json = ic.to_dict()

# Grab the userid from AWS so we can place outputs in the Scratch Bucket
# under the users prefix from our dask workers.
s3 = boto3.client('s3')
scratch_bucket = "adias-prod-user-scratch"
userid = boto3.client('sts').get_caller_identity()['UserId']

# Set this prefix to place all the data from this run under it in s3
s3_outputs_prefix = 'tasmania'


profiler = cProfile.Profile()
profiler.enable()

for item in items_json['features']:
    result = download_and_process(item, userid, scratch_bucket, s3_outputs_prefix)
    print(result)

profiler.disable()
stats = pstats.Stats(profiler)
stats.dump_stats(f'./{current_time}_stats_file.dat')
stats.strip_dirs()
stats.print_stats(15).sort_stats('cumtime')








# scene = "LC08_L1TP_114079_20220130_20220204_02_T1"
# # Each scene will have all of its assets download into a sub-directory with it's id.
# inputs_path = Path(f'/tmp/sample_data/{scene}/inputs')
# outputs_path = Path(f"/tmp/sample_data/{scene}/outputs")
# inputs_path.mkdir(parents=True, exist_ok=True)

# settings = {
#     "inputfile": str(inputs_path),
#     "output": str(outputs_path),
#     # polygon=
#     # limit=-29.9,152.5,-29.2,154.0
#     "l2w_parameters": "Rrs_*",
#     "rgb_rhot": False,
#     "rgb_rhos": False,
#     # "map_l2w": False,
#     # "merge_zones": False,
#     # "dsf_residual_glint_correction": True,
#     "l2w_export_geotiff": True,
#     "l1r_delete_netcdf" : True,
#     "l2r_delete_netcdf": True,
#     "l2t_delete_netcdf": True,
#     "l2w_delete_netcdf": True,
#     # "dsf_interface_reflectance": False, # False is the default
#     "ancillary_data": False, # If you set this to True you must supply a username and password for EARTHDATA
#     # "EARTHDATA_u": "",
#     # "EARTHDATA_p": "",
#     # "verbosity": 5
# }

# profiler = cProfile.Profile()
# profiler.enable()

# ac.acolite.acolite_run(settings=settings)

# profiler.disable()
# stats = pstats.Stats(profiler)
# stats.dump_stats(f'./{current_time}_stats_file.dat')
# stats.sort_stats('cumtime')
# stats.strip_dirs()
# stats.print_stats(15)