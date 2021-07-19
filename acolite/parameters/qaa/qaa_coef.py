## def qaa coeff
## reads qaa coefficients
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-08
## modifications:
##                2018-07-18 (QV) changed acolite import name
##                2021-03-29 (QV) new generic acolite, added bools

def qaa_coef():
    import os,sys
    import acolite as ac
    nfile = ac.config['data_dir']+'/Shared/algorithms/QAA/qaa_settings.txt'

    data = {}
    with open(nfile, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            if line[0] in [';','!', '/']: continue
            line = line.strip()
            if len(line)==0: continue
            split = line.split('=')
            if len(split)==2:
                if ',' in split[1]:
                    data[split[0]]=[float(d) for d in split[1].split(',')]
                else:
                    if split[1] in ['True', 'true']:
                        data[split[0]]=True
                    elif split[1] in ['False', 'false']:
                        data[split[0]]=False
                    else:
                        data[split[0]]=split[1]
    if 'useconfig' in data:
        for tag in ['g','h','k','l','m']:
            data[tag] = data['{}_{}'.format(tag,data['useconfig'])]
    return(data)
