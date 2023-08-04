## def read_6sv_ph
## reads 6SV aerosol model phase function from ACOLITE distribution, computes bbr
## written by Quinten Vanhellemont, RBINS
## 2023-07-03
## modifications: 2023-07-03 (QV) reduced the tolerance for the integration to 1e-4
##

def read_6sv_ph(model, compute_bbr=True):
    import acolite as ac
    import numpy as np
    import os

    import scipy.integrate as integrate

    ## function for integration of bbr
    def interp_ph(x):
        return(0.5 * np.interp(x, cangrad, cdata[:, wi]) * np.sin(x))

    model_dir = '{}/Shared/6SV/'.format(ac.config['data_dir'])
    models = [os.path.splitext(m)[0] for m in os.listdir(model_dir)]

    if model.upper() in ['1', 'C', 'MOD1', 'continental', 'continental_ph_6sv']:
        model_name = 'continental_ph_6sv'
        model_char = 'C'
        model_idx = 1
    elif model.upper() in ['2', 'M', 'MOD2', 'maritime', 'maritime_ph_6sv']:
        model_name = 'maritime_ph_6sv'
        model_char = 'M'
        model_idx = 2
    elif model.upper() in ['3', 'U', 'MOD3', 'urban', 'urban_ph_6sv']:
        model_name = 'urban_ph_6sv'
        model_char = 'U'
        model_idx = 3
    else:
        print('Model {} not recognised.'.format(model))
        return

    cfile = '{}/Shared/6SV/{}.csv'.format(ac.config['data_dir'], model_name)
    if os.path.exists(cfile):
        cwave = np.loadtxt(cfile, delimiter=',', max_rows=1, dtype=str)
        cwave = np.asarray([float(c.strip('"')) for c in cwave[1:]])

        cang = np.loadtxt(cfile, delimiter=',', usecols=[0], dtype=str)
        cang = np.asarray([float(c.strip('"')) for c in cang[1:]])
        cang = np.flip(cang)
        cangrad = np.radians(cang)

        cdata = np.loadtxt(cfile, delimiter=',', skiprows=1, usecols=[i+1 for i in range(len(cwave))], dtype=float)
        cdata = np.flipud(cdata)

        ## return bbr
        if compute_bbr:
            bbr = np.zeros(len(cwave))
            for wi, w in enumerate(cwave):
                #I = integrate.quadrature(interp_ph, np.pi/2, np.pi, maxiter=200)
                I = integrate.quad(interp_ph, np.pi/2, np.pi, epsabs=1e-4)
                bbr[wi] = I[0]
            return(bbr)

        return(cdata)
