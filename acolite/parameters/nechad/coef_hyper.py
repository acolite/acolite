## def coef_nechad_hs
## reads Nechad hyperspectral data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-05-15
## modifications:
##                2018-07-18 (QV) changed acolite import name
##                2021-03-09 (QV) adapted for acg
##                2023-07-30 (QV) added year, added "S" as alias for SPM

def coef_hyper(par, year = None):
    import os,sys
    import acolite as ac
    import numpy as np

    if par.upper() in ['S', 'SPM', 'TSM']:
        if year is None: year = 2010
        file = ac.config['data_dir']+'/Shared/algorithms/Nechad//SPM_N{}_Published.txt'.format(year)

    if par.upper() in ['T', 'TUR', 'TURBIDITY']:
        if year is None: year = 2009
        file = ac.config['data_dir']+'/Shared/algorithms/Nechad//Turbidity_N{}_Published.txt'.format(year)

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
