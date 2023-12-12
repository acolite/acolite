import acolite as ac
import keyring

InFile = '/D/Data/TEST/TEST_WISE/MC-10/190820_MC-10A-WI-1x1x1_v02-L1G.pix'
OutPath = '/D/Data/TEST/TEST_WISE/MC-10/L2_ACOLITE'

passwd = keyring.get_password("EarthData", "raphidoc")

Settings = {"inputfile":InFile,
            "output":OutPath,
            "polygon_limit":True,
            'aerosol_correction':'dark_spectrum',
            "dsf_residual_glint_correction": True,
            'dsf_write_aot_550':True,
            'dsf_write_tiled_parameters':False,
            'glint_write_rhog_ref':False,
            "l2w_parameters":['rhow_*'],
            #"output_xy":True,
            "output_geolocation":True,
            "output_xy":False,
            "reproject_outputs":'L2w',
            "output_projection":False,
            "reproject_before_ac":False,
            "output_projection_epsg":4326,
            'netcdf_projection':True,
            "ancillary_data":True,
            "EARTHDATA_u":"raphidoc",
            "EARTHDATA_p":passwd
            }

ac.acolite.acolite_run(settings=Settings)
