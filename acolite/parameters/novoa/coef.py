## def coef
## reads Novoa switching data
## written by Quinten Vanhellemont, RBINS
## 2023-07-30
## modifications:

def coef(config='defaults'):
    import os,sys
    import acolite as ac
    nfile = ac.config['data_dir']+'/Shared/algorithms/Novoa/{}.txt'.format(config)
    dct = {}
    with open(nfile, 'r') as f:
        for line in f.readlines():
            if line[0] in [';','!', '/', '#']: continue
            line = line.strip()
            split = line.split('=')
            if len(split) != 2: continue
            split[1] = split[1].split(',')
            try:
                dct[split[0]] = [float(s) for s in split[1]]
            except:
                dct[split[0]] = split[1]
            if len(dct[split[0]]) == 1: dct[split[0]]=dct[split[0]][0]
            if dct[split[0]] in ['True','true']: dct[split[0]]=True
            if dct[split[0]] in ['False','false']: dct[split[0]]=False
            if dct[split[0]] in ['None','none']: dct[split[0]]=None
    return(dct)
