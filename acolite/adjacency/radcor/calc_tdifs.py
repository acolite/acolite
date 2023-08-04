## def calc_tdifs
## computes diffuse transmittances
## written by Quinten Vanhellemont, RBINS
## 2023-06-14
## modifications:
##

def calc_tdifs(tau_ray, tau_aer, tau_gas, bbr = 0.0585, w0 = 1, theta_v = 0):
    import numpy as np
    B = bbr / w0
    cosv = np.cos(theta_v)

    tudir = np.exp(-(tau_ray + tau_aer + tau_gas) / cosv)
    tutot = np.exp(-(0.5 * tau_ray + B * tau_aer + tau_gas) / cosv)

    tuaer = tutot * np.exp(-0.5 * tau_ray / cosv) - tudir
    turay = tutot - tudir - tuaer

    return({'turay':turay, 'tuaer': tuaer, 'tudir': tudir, 'tutot': tutot})
