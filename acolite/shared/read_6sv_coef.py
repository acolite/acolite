## def read_6sv_coef
## reads 6SV aerosol model coefficients from ACOLITE distribution
## written by Quinten Vanhellemont, RBINS
## 2023-06-13
## modifications: 2023-06-13 (QV) added sensor keyword
##                2023-08-04 (QV) removed sensor keyword - perform convolution on the fly

def read_6sv_coef(model):
    import acolite as ac
    import numpy as np
    import os

    model_dir = '{}/Shared/6SV/'.format(ac.config['data_dir'])
    models = [os.path.splitext(m)[0] for m in os.listdir(model_dir)]

    if model.upper() in ['1', 'C', 'MOD1', 'continental', 'continental_coef_6sv']:
        model_name = 'continental_coef_6sv'
        model_char = 'C'
        model_idx = 1
    elif model.upper() in ['2', 'M', 'MOD2', 'maritime', 'maritime_coef_6sv']:
        model_name = 'maritime_coef_6sv'
        model_char = 'M'
        model_idx = 2
    elif model.upper() in ['3', 'U', 'MOD3', 'urban', 'urban_coef_6sv']:
        model_name = 'urban_coef_6sv'
        model_char = 'U'
        model_idx = 3
    else:
        print('Model {} not recognised.'.format(model))
        return

    cfile = '{}/Shared/6SV/{}.csv'.format(ac.config['data_dir'], model_name)
    if os.path.exists(cfile):
        cheader = np.loadtxt(cfile, delimiter=',', max_rows=1, dtype=str)
        cheader = np.asarray([c.strip('"') for c in cheader])
        cdata = np.loadtxt(cfile, delimiter=',', skiprows=1, dtype=float)
        cdata = {h: cdata[:, ih] for ih, h in enumerate(cheader)}
        return(cdata)
