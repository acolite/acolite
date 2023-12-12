## def l1_convert
## converts VIIRS data to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2023-03-30
## modifications: 2023-04-18 (QV) check if outside limits
##                2023-04-19 (QV) added quality_flags check
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call

def l1_convert(inputfile, output = None, settings = {}, verbosity = 0):
    import h5py
    import numpy as np
    import scipy.ndimage
    import dateutil.parser, time
    import glob, os
    import acolite as ac



    return(ofiles, setu)
