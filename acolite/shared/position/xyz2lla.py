## def xyz2lla
## compute lat, lon, alt from x,y,z coordinates by default using WGS84 axes (a & b)
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-20
## modifications:

def xyz2lla(x, y, z, a = 6378137.0, b = 6356752.314, degrees = True):
    import numpy as np

    r = np.sqrt(x**2 + y**2 + z**2)
    e2 = (a**2 - b**2) / a**2
    p = np.sqrt(x**2 + y**2)

    ## compute lat/lon
    lat = np.degrees(np.arctan(z / ((1- e2) * p)))
    lon = np.degrees(np.arctan(y/x))

    ## compute altitude
    v = a / np.sqrt(1 - e2 * (np.sin(np.radians(lat)))**2)
    alt = (p / np.cos(np.radians(lat))) - v

    return(lon, lat, alt)
