## def f0_get
## reads f0 data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-24
## modifications:
##                2017-11-28 (QV) moved PP data directory
##                2018-07-18 (QV) changed acolite import name
##                2018-09-19 (QV) added encoding
##                2021-02-05 (QV) changes for acolite-gen
##                2022-09-28 (QV) added support for bzipped files
##                2022-11-16 (QV) added printout when file is not found

def f0_get(f0_file=None, f0_dataset='Thuillier2003'):
    import numpy as np
    import os,sys,bz2
    import acolite as ac

    if f0_file is None:
        for ext in ['txt', 'txt.bz2']:
            f0_file = '{}/data/Solar/{}.{}'.format(ac.path, f0_dataset, ext)
            if os.path.exists(f0_file):
                break

    f0data=[]
    f0wave=[]
    lines = []
    if os.path.exists(f0_file):
        p, e = os.path.splitext(f0_file)
        if e == '.txt':
            with open(f0_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
        elif e == '.bz2':
            with bz2.open(f0_file, 'rb') as f:
                lines = [l.decode('utf-8') for l in f.readlines()]

        for line in lines:
                line = line.strip()
                if line[0] in ['#', '!', '/']: continue
                split = line.split(' ')
                if len(split) != 2: continue
                f0data.append(float(split[1]))
                f0wave.append(float(split[0]))
        f0={"wave":np.asarray(f0wave), "data":np.asarray(f0data)}
        return(f0)
    else:
        print('F0 reference {} not found in {}/data/Solar/'.format(f0_dataset, ac.path))
        return()
