## def parameter_cf_attributes
## gets cf attributes for parameters
##
## Written by Quinten Vanhellemont 2021-06-04
## Last modifications:

def parameter_cf_attributes():
    import acolite as ac
    import numpy as np
    param = {}
    header = None
    with open(ac.config['parameter_cf_attributes'], 'r', encoding="utf-8") as f:
        for line in f.readlines():
            line = line.strip()
            if len(line) == 0: continue
            if line[0] in ['#',';']: continue
            split = [i.strip() for i in line.split(',')]
            if header is None:
                header = split
            else:
                if header is None: continue
                par = split[0]
                val = {h:split[i] for i,h in enumerate(header)}
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
