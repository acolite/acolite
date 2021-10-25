## import sky reflectance luts for wind speed range
## QV 2021-01-31
## last updates: 2021-02-01 (QV) added sensor support
##               2021-02-24 (QV) renamed from rsky_read_luts
##               2021-03-01 (QV) simplified, added wind in lut
##               2021-10-24 (QV) added get_remote as keyword

def import_rsky_luts(models=[1,2], lutbase='ACOLITE-RSKY-202102-82W', sensor=None,
                    override=False, get_remote = True):
    import os
    import numpy as np
    import scipy.interpolate
    from netCDF4 import Dataset
    import acolite as ac

    rskyd = {}

    for mod in models:
        lut = None
        ret = ac.aerlut.import_rsky_lut(mod, lutbase=lutbase, sensor=sensor, get_remote = get_remote)
        rskyd[mod] = {'lut':ret[0], 'meta':ret[1], 'dims':ret[2], 'rgi':ret[3]}

        if sensor is None: ## generic model
            rskyd[mod]['dim'] = [rskyd[mod]['meta']['wave'], rskyd[mod]['meta']['azi'], rskyd[mod]['meta']['thv'],
                                 rskyd[mod]['meta']['ths'], rskyd[mod]['meta']['wind'], rskyd[mod]['meta']['tau']]
            rskyd[mod]['rgi'] = scipy.interpolate.RegularGridInterpolator(rskyd[mod]['dim'],
                                                                          rskyd[mod]['lut'],
                                                                          bounds_error=False,
                                                                          fill_value=np.nan)
        else: ## sensor specific model
            rskyd[mod]['dim'] = [rskyd[mod]['meta']['azi'], rskyd[mod]['meta']['thv'],
                                 rskyd[mod]['meta']['ths'], rskyd[mod]['meta']['wind'], rskyd[mod]['meta']['tau']]

            rskyd[mod]['rgi'] = {}
            for b in rskyd[mod]['lut']:
                rskyd[mod]['rgi'][b] = scipy.interpolate.RegularGridInterpolator(rskyd[mod]['dim'],
                                                                                 rskyd[mod]['lut'][b],
                                                                                 bounds_error=False,
                                                                                 fill_value=np.nan)
    return(rskyd)
