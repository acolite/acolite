## def coef_nechad_hs
## reads Nechad hyperspectral data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-05-15
## modifications:
##                2018-07-18 (QV) changed acolite import name
##                2021-03-09 (QV) adapted for acg

def coef_hyper(par):
    import os,sys
    import acolite as ac
    import numpy as np

    if par.upper() in ['SPM', 'TSM']:
        file = ac.config['data_dir']+'/Shared/algorithms/Nechad//SPM_N2010_Published.txt'

    if par.upper() in ['T', 'TUR', 'TURBIDITY']:
        file = ac.config['data_dir']+'/Shared/algorithms/Nechad//Turbidity_N2009_Published.txt'

    keys = ['wave','A','B','Rsq','C']
    data = {k:[] for k in keys}
    with open(file,'r') as f:
        for line in f.readlines():
            if line[0] in ['#',';', '%']: continue
            sp = line.strip().split(',')
            if len(sp) != 5: continue
            for i,k in enumerate(keys):
                data[k].append(float(sp[i]))
    for k in data: data[k] = np.asarray(data[k])

    return(data)
