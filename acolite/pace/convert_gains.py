## def convert_gains
## function to convert PACE OCI NetCDF gains file from OCSSW
## from 286 to 291 bands length
##
## last three bands for the blue detector are set to 1,0
## values are replicated for SWIR high/low gain bands
##
## written by Quinten Vanhellemont, RBINS
## 2026-06-02
## modifications:
##

def convert_gains(gains_nc):
    import xarray as xr
    import numpy as np
    import acolite as ac

    ## read gains
    with xr.open_dataset(gains_nc) as ds:
        gain_wave = ds['wavelength'].values
        gain = ds['gain'].values
        offset = ds['offset'].values

    ## reformat to all TOA bands
    ## create empty gains data
    oci_l1b_bandpass = ac.shared.obpg_bandpass('oci_l1b')
    toa_gains = np.zeros(len(oci_l1b_bandpass[:, 2]))
    toa_offsets = np.zeros(len(oci_l1b_bandpass[:, 2]))

    ## copy first bands
    toa_gains[0:116] = gain[0:116]
    toa_gains[116:119] = 1.0 ## set last three blue bands to one
    toa_gains[119:282] = gain[116:279]

    toa_offsets[0:116] = offset[0:116]
    toa_offsets[116:119] = 0.0 ## set last three blue bands to zero
    toa_offsets[119:282] = offset[116:279]

    ## copy swir bands
    swir_bands_in = [279, 280, 281, 281, 282, 283, 283, 284, 285]
    swir_bands_out = range(282, 291)
    for i, o in enumerate(swir_bands_out):
        toa_gains[o] = gain[swir_bands_in[i]]
        toa_offsets[o] = offset[swir_bands_in[i]]

    return(toa_gains, toa_offsets)
