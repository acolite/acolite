## def gas_transmittance
## gets hyperspectral gas transmittances
## written by Quinten Vanhellemont, RBINS
## 2019-04-02
## modifications:

def gas_transmittance(sza, vza, waves = None, uoz = 0.3, uwv = 1.5,
                      rsr = None, sensor = None):
    import acolite as ac
    import numpy as np

    ko3 = ac.ac.ko3_read()
    if waves is None: waves=ko3['wave']
    else: waves = [float(w/1000) for w in waves]

    ## compute water transmittance
    wv_wv_hs, tt_wv_hs = ac.ac.wvlut_interp(sza, vza, uwv=uwv)
    tt_wv = np.interp(waves, wv_wv_hs, tt_wv_hs)

    ## compute oxygen transmittance
    wv_o2_hs, tt_o2_hs = ac.ac.o2lut_interp(sza, vza)
    tt_o2 = np.interp(waves, wv_o2_hs, tt_o2_hs)

    ## cosine of sun and sensor zenith angles
    mu0 = np.cos(sza*(np.pi/180))
    muv = np.cos(vza*(np.pi/180))

    ## compute ozone transmittance
    koz = np.interp(waves, ko3['wave'], ko3['data'])
    tau_oz = koz * uoz
    t0_ozone = np.exp(-1.*(tau_oz) / mu0)
    tv_ozone = np.exp(-1.*(tau_oz) / muv)
    tt_o3= t0_ozone * tv_ozone

    ## prepare gas transmittance dict
    d =  {'wave':[1000*w for w in waves],
            'tt_h2o':tt_wv,
            'tt_o3':tt_o3,
            'tt_o2':tt_o2, 'tt_gas':tt_wv*tt_o3*tt_o2}

    ## resample if sensor is requested
    if sensor is not None:
        # find RSR
        rsr_file = ac.config['data_dir']+'/RSR/{}.txt'.format(sensor)
        rsr,bands = ac.shared.rsr_read(file=rsr_file)

    ## resample individual datasets
    if rsr is not None:
        for k in d: d[k] = ac.shared.rsr_convolute_dict(waves, d[k], rsr)

    return(d)
