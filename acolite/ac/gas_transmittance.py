## def gas_transmittance
## gets hyperspectral gas transmittances
## written by Quinten Vanhellemont, RBINS
## 2019-04-02
## modifications: 2022-11-17 (QV) added 2D interpolation for a given sensor, added new tto3_interp function

def gas_transmittance(sza, vza, pressure = 1013, waves = None, uoz = 0.3, uwv = 1.5,
                      gases = ['h2o', 'o3', 'o2', 'co2', 'n2o', 'ch4'], lutconfig = '202106F',
                      rsr = None, sensor = None):
    import acolite as ac
    import numpy as np

    ## hyperspectral
    if sensor is None:
        ## compute ozone transmittance
        wv_o3_hs, tt_o3_hs = ac.ac.tto3_interp(sza, vza, uoz=uoz)

        ## compute co2, o2, n2o, ch4 transmittance
        tg_hyp = ac.ac.gaslut_interp(sza, vza, pressure=pressure, lutconfig=lutconfig)

        ## compute water transmittance
        wv_wv_hs, tt_wv_hs = ac.ac.wvlut_interp(sza, vza, uwv=uwv)

        ## interpolate to same wavelengths
        if waves is None:
            waves = wv_o3_hs * 1.0
        else:
            waves = np.asarray([float(w/1000) for w in waves])
        tt_o3 = np.interp(waves, wv_o3_hs, tt_o3_hs)
        tg_wv_hs = tg_hyp['wave'] * 1.0
        for k in tg_hyp: tg_hyp[k] = np.interp(waves, tg_wv_hs, tg_hyp[k])
        tt_wv = np.interp(waves, wv_wv_hs, tt_wv_hs)

        ## prepare gas transmittance dict
        d =  {'wave': [1000*w for w in waves],
                'tt_h2o': tt_wv,
                'tt_o3': tt_o3,
                'tt_o2': tg_hyp['ttoxyg'], 'tt_co2': tg_hyp['ttdica'],
                'tt_n2o': tg_hyp['ttniox'], 'tt_ch4': tg_hyp['ttmeth']}

        ## total gas transmittance
        d['tt_gas'] = np.ones(len(d['wave']))
        for g in gases: d['tt_gas'] *= d['tt_{}'.format(g)]

        ## resample to rsr
        if rsr is not None:
            for k in d: d[k] = ac.shared.rsr_convolute_dict(waves, d[k], rsr)
    else:
        ## compute co2, o2, n2o, ch4 transmittance
        tg_ret = ac.ac.gaslut_interp(sza, vza, pressure=pressure, lutconfig=lutconfig, sensor=sensor)

        ## compute water transmittance
        wv_ret = ac.ac.wvlut_interp(sza, vza, uwv=uwv, sensor=sensor)

        ## compute ozone transmittance
        tt_o3 = ac.ac.tto3_interp(sza, vza, uoz=uoz, sensor=sensor)

        ## prepare gas transmittance dict
        d =  {'tt_h2o': wv_ret, 'tt_o3': tt_o3,
              'tt_o2': tg_ret['ttoxyg'], 'tt_co2': tg_ret['ttdica'],
              'tt_n2o': tg_ret['ttniox'], 'tt_ch4': tg_ret['ttmeth']}

        ## total gas transmittance
        d['tt_gas'] = {}
        for b in d['tt_h2o']:
            d['tt_gas'][b] = np.ones(d['tt_h2o'][b].shape)
            for g in gases: d['tt_gas'][b] *= d['tt_{}'.format(g)][b]
    return(d)
