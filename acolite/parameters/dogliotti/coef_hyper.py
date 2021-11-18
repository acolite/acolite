## def hyper
## reads new Dogliotti hyperspectral calibration data
## written by Quinten Vanhellemont, RBINS
## 2021-11-18
## modifications:

def coef_hyper():
    import os,sys
    import acolite as ac
    import numpy as np

    dfile = ac.config['data_dir'] + '/Shared/algorithms/Dogliotti/dogliotti_turbidity_hyperspectral.csv'

    ## read header
    with open(dfile, 'r', encoding='utf-8') as f:
        header = f.readline()
        header = header[1:].strip().split(',')
    ## read data
    tmp = np.loadtxt(dfile, delimiter=',')
    tmp.shape

    dogliotti = {}
    for ih, h in enumerate(header):
        dogliotti[h] = tmp[:, ih]
    tmp = None

    return(dogliotti)
