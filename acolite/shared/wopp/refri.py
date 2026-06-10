## def wopp.refri
## load WOPP refractive index of water for 27°C and 0 PSU
##
## Röttgers, R., R. Doerffer, D. McKee, and W. Schönfeld.
## "The water optical properties processor (WOPP): pure water spectral
## absorption, scattering and real part of refractive index model.""
##
## WOPP is available from https://calvalportal.ceos.org/tools
## txt files at data/Shared/WOPP/ were converted from .dat files
##
## written by Quinten Vanhellemont, RBINS
## 2018-07-30
## modifications: 2021-03-30 (QV) new acolite
##                2026-06-10 (QV) moved to shared.wopp.refri

def refri():
    import acolite as ac
    import numpy as np

    file = '{}/{}'.format(ac.config['data_dir'], 'Shared/WOPP/computed_refri_T27_S0.txt')

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
