## def ko3_read
## reads ko3 data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-04-14
## modifications:
##                2017-11-28 (QV) moved PP data directory
##                2018-07-18 (QV) changed acolite import name
##                2021-02-24 (QV) renamed from ko3_get

def ko3_read(ko3file=None):
    import os,sys
    import numpy as np
    import acolite as ac

    if ko3file is None: ko3file = ac.config['data_dir']+'/Shared/k_o3_anderson.txt'
    ko3data=[]
    ko3wave=[]
    with open(ko3file, 'r') as f:
        for line in f:
            if line[0] == '!': continue
            if line[0] == '/': continue
            split = line.split(' ')
            if len(split) != 2: continue
            ko3data.append(float(split[1]))
            ko3wave.append(float(split[0])/1000.)
    ko3={"wave":np.asarray(ko3wave), "data":np.asarray(ko3data)}
    return(ko3)
