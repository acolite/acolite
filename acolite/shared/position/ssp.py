## def ssp
## computes sub satelite point based on x,y,z
##
## written by Quinten Vanhellemont, RBINS
## based on notes from Satellite Times by Dr. T.S. Kelso (https://celestrak.com/columns/)
## 2025-05-21
## modifications: 2025-05-22 (QV) split from eci_geometry

def ssp(sat_x, sat_y, sat_z, time, first_estimate = False, Re = 6378.135, a = 6378.137, f = 1/298.257223563, ssp_tolerance = 0.001):
    import acolite as ac
    import numpy as np

    ## get theta_g
    theta_g = ac.shared.position.gmst(time)

    ## first estimate of SSP
    slat = np.arctan(sat_z / np.sqrt(sat_x**2 + sat_y**2))
    slon = np.pi + np.arctan(sat_y / sat_x) - theta_g
    salt = np.sqrt(sat_x**2 + sat_y**2 + sat_z**2) - Re
    if first_estimate: return(np.degrees(slon), np.degrees(slat), salt)

    ## iterate for exact SSP latitude and height
    esq = 2 * f - f ** 2
    R = np.sqrt(sat_x**2 + sat_y**2)
    diff = np.inf
    ## only for one point!
    while diff > ssp_tolerance:
        slat0 = 1 * slat
        C = 1 / np.sqrt(1 - esq * np.sin(slat) ** 2)
        slat = np.arctan((sat_z + a * C * esq * np.sin(slat))/R)
        diff = np.abs(slat - slat0)
    ## satellite height
    salt = R / np.cos(slat) - a * C

    return(np.degrees(slon), np.degrees(slat), salt)
