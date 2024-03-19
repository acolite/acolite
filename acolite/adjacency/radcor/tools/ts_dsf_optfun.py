## def ts_dsf_optfun
## optimisation function for TS-DSF AOT estimation
##
## written by Alexandre Castagna, UGent (R version)
##            Quinten Vanhellemont, RBINS (Python version)
##
## for the RAdCor project 2023-2024
##
## modifications: 2024-03-13 (QV) added as function

def ts_dsf_optfun(x, a_min, b_min, lut_rho_a, lut_T_u_dif_r):
    import acolite as ac
    import numpy as np

    est = (a_min - x) / (b_min - x)
    exp = np.interp(x, lut_rho_a, lut_T_u_dif_r)
    cost = (exp - est)**2
    if ac.settings['run']['radcor_development']:
        print('        Test rho_a: {:.6f}  Est. ratio: {:.5f}  Exp. ratio: {:.5f}  Cost value: {:.3E}'.format(x, est, exp, cost))
    return(cost)
