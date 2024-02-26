## def irradiance
## compute surface level irradiance for given sun zenith angle, aot and aerosol model
##
## QV 2023-05-31
## modifications: 2024-02-26 (QV) added downward gas transmittance, pressure, multiplication with cosine sun zenith angle

def irradiance(aot, sza, mod = 2, date = None, apply_cosine = True, settings = {}):
    import acolite as ac
    import numpy as np
    import dateutil

    ## parse settings
    setu = ac.acolite.settings.parse(None, settings=settings, merge=True)

    ## get defaults
    pressure = setu['pressure']
    uoz = setu['uoz_default']
    uwv = setu['uwv_default']

    ## read f0
    f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
    wave = f0['wave']/1000
    raa=0
    vza=0

    ## sun earth distance correction
    if date is None:
        dist = 1
    else:
        dtime = dateutil.parser.parse(date)
        doy = dtime.strftime('%j')
        dist = ac.shared.distance_se(doy)
    dsq = dist**2
    f0['data'] *= dsq

    ## get gas transmittance
    ## old version: two-way which is not correct for surface irradiance
    #tg_dict = ac.ac.gas_transmittance(sza, vza, uoz=uoz, uwv=uwv)

    ## new version with downward path only
    ## compute co2, o2, n2o, ch4 transmittance
    tg_hyp = ac.ac.gaslut_interp(sza, vza, pressure=pressure, lutconfig='202402F',
                            pars = ['dtdica','dtoxyg','dtniox','dtmeth'])

    ## compute ozone transmittance
    tg_o3 = ac.ac.tto3_interp(sza, vza, uoz=uoz, total=False)

    ## interpolate to same grid
    tg_hyp['dto3'] = np.interp(tg_hyp['wave'], tg_o3['wave'], tg_o3['dt_o3'])

    ## compute water transmittance
    wv_wv_hs, dt_wv_hs = ac.ac.wvlut_interp(sza, vza, uwv=uwv, pressure=pressure, par = 'dtwava', config = '202402CP')
    tg_hyp['dth2o'] = dt_wv_hs

    ## total gas transmittance
    gases = [k for k in tg_hyp if k[0:2] == 'dt']
    tg_hyp['dtgas'] = np.ones(len(tg_hyp['wave']))
    for g in gases: tg_hyp['dtgas'] *= tg_hyp[g]

    ## interpolate to f0 grid
    tgas = np.interp(f0['wave']/1000, tg_hyp['wave'], tg_hyp['dtgas'])

    ## load LUT
    lutid = [lut for lut in setu['luts'] if 'MOD{}'.format(mod) in lut][0]
    lutdw = ac.aerlut.import_luts(add_rsky=False, par='romix', sensor=None, base_luts=[lutid])

    ## compute downward total transmittance
    res = lutdw[lutid]['rgi']((pressure, lutdw[lutid]['ipd']['dtott'], wave, raa, vza, sza, aot))

    ## compute Ed
    ed = res * f0['data'] * tgas

    if apply_cosine:
        ed *= np.cos(np.radians(sza))

    dct = {'wave': wave*1000, 'ed': ed, 'f0': f0['data'], 'dtott': res, 'tgas': tgas}
    return(dct)
