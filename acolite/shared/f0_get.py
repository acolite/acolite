## def f0_get
## reads f0 data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-24
## modifications:
##                2017-11-28 (QV) moved PP data directory
##                2018-07-18 (QV) changed acolite import name
##                2018-09-19 (QV) added encoding
##                2021-02-05 (QV) changes for acolite-gen

def f0_get(f0_file=None, f0_dataset='Thuillier2003'):
    import numpy as np

    if f0_file is None:
        import os,sys
        import acolite as ac
        f0_file = '{}/data/Solar/{}.txt'.format(ac.path, f0_dataset)

    f0data=[]
    f0wave=[]
    with open(f0_file, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            line = line.strip()
            if line[0] in ['#', '!', '/']: continue
            split = line.split(' ')
            if len(split) != 2: continue
            f0data.append(float(split[1]))
            f0wave.append(float(split[0]))
    f0={"wave":np.asarray(f0wave), "data":np.asarray(f0data)}
    return f0
