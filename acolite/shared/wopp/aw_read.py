## def wopp.aw_read
## function to read WOPP aw data
##
## Röttgers, R., R. Doerffer, D. McKee, and W. Schönfeld.
## "The water optical properties processor (WOPP): pure water spectral
## absorption, scattering and real part of refractive index model.""
##
## WOPP is available from https://calvalportal.ceos.org/tools
## txt files at data/Shared/WOPP/ were converted from .dat files
##
## written by Quinten Vanhellemont, RBINS
## 2021-09-22
## modifications:

def aw_read(coef_file = None, version = 1):
    import acolite as ac
    import numpy as np

    ## purewater file from WOPP
    if coef_file is None:
        coef_file = '{}/{}'.format(ac.config['data_dir'], 'Shared/WOPP/purewater_abs_coefficients_v{}.txt'.format(version))

    data = {'wave':[], 'a':[], 'PsiS':[], 'PsiT':[], 'sigma_a':[], 'sigma_PsiS':[], 'sigma_PsiT':[]}
    with open(coef_file, 'r', encoding = "utf-8") as f:
        for line in f.readlines():
            if line[0] in ['%', '#', ';']: continue
            line = line.strip()
            s = line.split()
            if len(s) == 7:
                data['wave'].append(float(s[0]))
                data['a'].append(float(s[1]))
                data['PsiS'].append(float(s[2]))
                data['PsiT'].append(float(s[3]))
                data['sigma_a'].append(float(s[4]))
                data['sigma_PsiS'].append(float(s[5]))
                data['sigma_PsiT'].append(float(s[6]))
    for k in data: data[k] = np.asarray(data[k])
    return(data)
