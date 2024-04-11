import acolite as ac
import keyring

#InFile = '/D/Data/TEST/TEST_WISE/MC-10/190820_MC-10A-WI-1x1x1_v02-L1G.pix'
# After L1R conversion
InFile = '/D/Data/TEST/TEST_WISE/MC-10/L2_ACOLITE/wise_2606_190820_MC-10A-WI-1x1x1_v02_2019_08_20_15_11_38_L1R.nc'
OutPath = '/D/Data/TEST/TEST_WISE/MC-10/L2_ACOLITE'

passwd = keyring.get_password("EarthData", "raphidoc")

Settings = {"inputfile":InFile,
            "output":OutPath,
            "blackfill_skip":True,
            "blackfill_wave":802,
            'aerosol_correction':'dark_spectrum',
            'dsf_write_aot_550':True,
            'dsf_write_tiled_parameters':False,
            "dsf_residual_glint_correction": False,
            'glint_write_rhog_ref':False,
            "l2w_parameters":['rhow_*'],
            "l2w_mask_high_toa":False,
            "output_geolocation":True,
            "output_xy":True,
            "output_projection":False,
            'netcdf_projection':True,
            "ancillary_data":True,
            "EARTHDATA_u":"raphidoc",
            "EARTHDATA_p":passwd,
            "output_lt":False,
            "l2w_mask_negative_rhow":False
            }

ac.acolite.acolite_run(settings=Settings)
