## def o3_interp
## gets ozone transmittances
## written by Quinten Vanhellemont, RBINS
## 2022-11-17 (split off from gas_transmittance)
## modifications: 2024-02-26 (QV) added total keyword

def tto3_interp(sza, vza, uoz = 0.3, sensor = None, total = True):

    import acolite as ac
    import numpy as np

    ## import ozone pars
    ko3 = ac.ac.ko3_read()

    ## cosine of sun and sensor zenith angles
    mu0 = np.cos(sza*(np.pi/180))
    muv = np.cos(vza*(np.pi/180))

    ## compute ozone transmittance
    if sensor is None:
        tau_oz = ko3['data'] * uoz
        t0_ozone = np.exp(-1.*(tau_oz) / mu0)
        tv_ozone = np.exp(-1.*(tau_oz) / muv)
        tt_o3= t0_ozone * tv_ozone
        if total:
            return(ko3['wave'], tt_o3)
        else:
            return({'wave':ko3['wave'], 'tt_o3': tt_o3, 'dt_o3': t0_ozone, 'ut_o3': tv_ozone})
    else:
        ## find RSR
        rsrd = ac.shared.rsr_dict(sensor=sensor)
        if sensor in rsrd:
            rsr, rsr_bands = rsrd[sensor]['rsr'], rsrd[sensor]['rsr_bands']
        else:
            print('Sensor {} RSR not found'.format(sensor))
            sys.exit(1)

        ## compute ozone transmittance
        koz = ac.shared.rsr_convolute_dict(ko3['wave'], ko3['data'], rsr)
        tt_o3 = {}
        dt_o3 = {}
        ut_o3 = {}
        for b in koz:
            tau_oz = koz[b] * uoz
            t0_ozone = np.exp(-1.*(tau_oz) / mu0)
            tv_ozone = np.exp(-1.*(tau_oz) / muv)
            tt_o3[b] = t0_ozone * tv_ozone
            dt_o3[b] = t0_ozone
            ut_o3[b] = tv_ozone

        if total:
            return(tt_o3)
        else:
            return({'tt_03': tt_o3, 'dt_o3': dt_o3, 'ut_o3': ut_o3})
