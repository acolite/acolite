## def tss
## function to compute Jiang et al. 2021 (OLCI) or 2023 (MSI) TSS for Rrs in a given band
##
## written by Quinten Vanhellemont, RBINS
## 2025-06-30
## modifications:

def tss(Rrs, tss_cfg, band_wave, sub = None):
    import numpy as np

    ## function for subsurface rrs
    def _rrs(Rrs, sub = None):
        if sub is None:
            return(Rrs/(0.52+1.7*Rrs))
        else:
            return(Rrs[sub]/(0.52+1.7*Rrs[sub]))

    ## get band index for aw and bbw
    w_idx = np.argsort(np.abs(np.asarray(tss_cfg['Rrs_wave']) - band_wave))[0]
    ## get band index for tss
    t_idx = np.argsort(np.abs(np.asarray(tss_cfg['TSS_wave']) - band_wave))[0]

    ## QAA computation
    ## subsurface rrs
    rrs = _rrs(Rrs[band_wave], sub = sub)
    ## u
    u = (-0.089+((0.089**2)+4*0.125*rrs)**0.5)/(2*0.125)

    ## bbp
    if band_wave in [560.0, 665.0]:
        if band_wave == 560.0:
            rrs443 = _rrs(Rrs[443.0], sub = sub)
            rrs490 = _rrs(Rrs[490.0], sub = sub)
            rrs665 = _rrs(Rrs[665.0], sub = sub)
            x = np.log10((rrs443+rrs490)/(rrs+5*rrs665*rrs665/rrs490))
            del rrs443, rrs490, rrs665
            a = tss_cfg['aw'][w_idx]+10**(-1.146-1.366*x-0.469*(x**2))
        elif band_wave == 665.0:
            if sub is None:
                a = tss_cfg['aw'][w_idx]+0.39*((Rrs[665.0]/(Rrs[443.0]+Rrs[490.0]))**1.14)
            else:
                a = tss_cfg['aw'][w_idx]+0.39*((Rrs[665.0][sub]/(Rrs[443.0][sub]+Rrs[490.0][sub]))**1.14)
        bbp = ((u*a)/(1-u))-tss_cfg['bbw'][w_idx]
        del a
    else:
        bbp = ((u*tss_cfg['aw'][w_idx])/(1-u))-tss_cfg['bbw'][w_idx]
    del u
    
    ## TSS computation
    tss = tss_cfg['TSS_scale'][t_idx] * bbp
    del bbp

    return(tss)
