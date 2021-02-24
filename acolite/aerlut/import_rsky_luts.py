## import sky reflectance luts for wind speed range
## QV 2021-01-31
## last updates: 2021-02-01 (QV) added sensor support
##               2021-02-24 (QV) renamed from rsky_read_luts

def import_rsky_luts(models=[1,2], lutbase='ACOLITE-RSKY-202101-75W-{}ms', winds=[1,2,5,10], sensor=None, override=False):
    import os
    import numpy as np
    import scipy.interpolate
    from netCDF4 import Dataset
    import acolite as ac

    rskyd = {}
    winds = np.asarray(winds)

    if sensor is None:
        for mod in models:
            lut = None
            for iw, wind in enumerate(winds):
                ret = ac.aerlut.import_rsky_lut(mod, lutbase=lutbase.format(wind))
                ret = {'lut':ret[0], 'meta':ret[1], 'dims':ret[2], 'rgi':ret[3]}
                tmp = ret['lut']
                if lut is None:
                    lut = tmp[:, :, :, :, None, :]
                else:
                    lut = np.insert(lut,(iw),tmp,axis=4)

            dim = [ret['meta']['wave'], ret['meta']['azi'], ret['meta']['thv'],
                   ret['meta']['ths'], winds, ret['meta']['tau']]

            rskyd[mod]={'lut':lut, 'meta':ret['meta'], 'dim':dim}
            rskyd[mod]['rgi'] = scipy.interpolate.RegularGridInterpolator(rskyd[mod]['dim'],
                                                                          rskyd[mod]['lut'],
                                                                          bounds_error=False,
                                                                          fill_value=np.nan)
    else:
        for mod in models:
            lut = None
            for iw, wind in enumerate(winds):
                ret = ac.aerlut.import_rsky_lut(mod, lutbase=lutbase.format(wind), sensor=sensor)
                ret = {'lut':ret[0], 'meta':ret[1], 'dims':ret[2], 'rgi':ret[3]}
                tmp = ret['lut']
                if lut is None:
                    lut = {b: tmp[b][:, :, :, None, :] for b in tmp}
                else:
                    for b in lut:
                        lut[b] = np.insert(lut[b],(iw),tmp[b],axis=3)

            dim = [ret['meta']['azi'], ret['meta']['thv'],
                   ret['meta']['ths'], winds, ret['meta']['tau']]

            rskyd[mod]={'lut':lut, 'meta':ret['meta'], 'dim':dim}
            rskyd[mod]['rgi'] = {}
            for b in rskyd[mod]['lut']:
                rskyd[mod]['rgi'][b] = scipy.interpolate.RegularGridInterpolator(rskyd[mod]['dim'],
                                                                          rskyd[mod]['lut'][b],
                                                                          bounds_error=False,
                                                                          fill_value=np.nan)
    return(rskyd)
