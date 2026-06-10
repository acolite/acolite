## def wopp.aw_ts
## function compute aw for given Temperature T1 and Salinity S1
## based on WOPP data
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

def aw_ts(T1, S1, T0 = 20, S0 = 0, aw = None):
    import acolite as ac
    if aw is None: aw = ac.shared.wopp.aw_read()

    ## copy wavelength to output dict
    data = {'wave': aw['wave']}
    data['awT'] = aw['a'] + (T1-T0) * aw['PsiT']
    data['awS'] = aw['a'] + (S1-S0) * aw['PsiS']

    data['awTS'] = aw['a'] + (T1-T0) * aw['PsiT'] + (S1-S0) * aw['PsiS']
    data['awTSe'] = aw['sigma_a'] + (T1-T0) * aw['sigma_PsiT'] + (S1-S0) * aw['sigma_PsiS']
    return(data)
