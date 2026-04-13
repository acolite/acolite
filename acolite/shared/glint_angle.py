## computes the angle between view vector and specular glint vector
## should be > 40 for low glint risk
##
## written by Quinten Vanhellemont, RBINS
## 2026-04-12
## modifications:

def glint_angle(sza, vza, raa, degrees_in = True, degrees_out = True):
    import numpy as np

    if degrees_in:
        theta_s = np.radians(sza)
        theta_v = np.radians(vza)
        delta_phi = np.radians(raa)
    else:
        theta_s = 1.0 * sza
        theta_v = 1.0 * vza
        delta_phi = 1.0 * raa

    ga = np.arccos(np.cos(theta_s) * np.cos(theta_v) + \
         np.sin(theta_s) * np.sin(theta_v) * np.cos(delta_phi))

    if degrees_out: ga = np.degrees(ga)

    return(ga)
