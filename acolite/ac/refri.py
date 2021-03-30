## load WOPP refractive index of water for 27°C and 0 PSU
## QV 2018-07-30
## modifications: 2021-03-30 (QV) new acolite

def refri():
    import acolite as ac
    import numpy as np

    file = '{}/{}'.format(ac.config['data_dir'], 'Shared/WOPP/computed_refri_T27_S0.dat')

    data = {'wave':[], 'n':[]}
    with open(file, 'r', encoding="utf-8") as f:
        for line in f.readlines():
            if line[0] in ['%', '#', ';']: continue
            line = line.strip()
            s = line.split()
            if len(s) == 2:
                data['wave'].append(float(s[0]))
                data['n'].append(float(s[1]))

    for k in data: data[k] = np.asarray(data[k])
    return(data)
