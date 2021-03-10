## def coef
## reads Dogliotti 2015 algorithm calibration data
## written by Quinten Vanhellemont, RBINS
## 2021-03-10
## modifications:

def coef(config='defaults'):
    import os,sys
    import acolite as ac
    nfile = ac.config['data_dir']+'/Shared/algorithms/Dogliotti/{}.txt'.format(config)
    dct = {}
    with open(nfile, 'r') as f:
        for line in f.readlines():
            if line[0] in [';','!', '/', '#']: continue
            line = line.strip()
            split = line.split('=')
            if len(split) != 2: continue
            try:
                dct[split[0]] = float(split[1])
            except:
                dct[split[0]] = split[1]
    return(dct)
