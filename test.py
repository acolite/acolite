import acolite as ac
import keyring

#InFile = '/D/Data/WISE/MC-10A/190820_MC-10A-WI-1x1x1_v02-L1G.dpix/190820_MC-10A-WI-1x1x1_v02-L1G.pix'
# After L1R conversion
InFile = '/D/Data/TEST/TEST_WISE/MC-10/L2_ACOLITE/wise_2606_190820_MC-10A-WI-1x1x1_v02_2019_08_20_15_11_38_L1R.nc'
OutPath = '/D/Data/TEST/TEST_WISE/MC-10/L2_ACOLITE'

# InFiles = ['/D/Data/TEST/merge/S2A_MSIL1C_20230824T152631_N0509_R068_T20UMA_20230824T203721.SAFE.zip',
#            '/D/Data/TEST/merge/S2A_MSIL1C_20230824T152631_N0509_R068_T20ULA_20230824T203721.SAFE']
# OutPath = '/D/Data/TEST/merge/test/'
#
# PolygonLimit = '/D/Data/TEST/merge/roi_acolite.geojson'

# PolygonLimit = os.path.join(r.Chapter2Path, "polygon_limit.geojson")
#
# with open(PolygonLimit, 'w') as outfile:
#     geojson.dump(bbox_json, outfile)


Settings = {"inputfile":InFile,
            "output":OutPath,
            #"merge_tiles":True,
            #"limit":[49.664140118419056, -64.4735034754853, 49.83986356624047, -64.08483698352762],
            #"polygon": PolygonLimit,
            #"polygon_limit": True,
            "blackfill_skip":True,
            #"blackfill_wave":802,
            'aerosol_correction':'dark_spectrum',
            'dsf_write_aot_550':True,
            'dsf_write_tiled_parameters':False,
            "dsf_residual_glint_correction": False,
            'glint_write_rhog_ref':False,
            "l2w_parameters":['rhow_*'],
            "l2w_mask_high_toa":False,
            "output_geolocation":True,
            "output_xy":True,
            "output_projection":True,
            'netcdf_projection':True,
            "ancillary_data":False,
            "output_lt":False,
            "l2w_mask_negative_rhow":False,
            "geometry_type":'grids'
            }

ac.acolite.acolite_run(settings=Settings)
