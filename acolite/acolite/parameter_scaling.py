## def pscale
## gets parameter scaling
##
## Written by Quinten Vanhellemont 2017-12-05
## Last modifications: 2018-04-17 (QV) added nan test
##                2018-07-18 (QV) changed acolite import name
##                2018-09-10 (QV) added encoding
##                2021-03-15 (QV) adapted for acg

def parameter_scaling():
    import acolite as ac
    import numpy as np
    param = {}
    header = None
    with open(ac.config['parameter_labels'], 'r', encoding="utf-8") as f:
        for line in f.readlines():
            line = line.strip()
            if len(line) == 0: continue
            if line[0] in ['#',';']: continue
            split = line.split('=')
            if split[0] == 'header':
                header = split[1].split(',')
            else:
                if header is None: continue
                tmp = [i.strip() for i in line.split(',')]
                par = tmp[0]
                val = {h:tmp[i] for i,h in enumerate(header)}
                for i in val:
                    if val[i] in ['False', 'false']: val[i] = False
                    if val[i] in ['True', 'true']: val[i] = True
                    if val[i] in ['None', 'none']: val[i] = None
                    if i in ['max', 'min']:
                        try:
                            val[i] = float(val[i])
                        except:
                            val[i] = np.nan
                param[par]=val

    return(param)
