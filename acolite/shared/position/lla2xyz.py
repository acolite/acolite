## def lla2xyz
## compute x,y,z coordinates from lat, lon, alt by default using WGS84 axes (a & b)
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-20
## modifications:

def lla2xyz(lon, lat, alt, a = 6378137.0, b = 6356752.314):
    import numpy as np

    lat = np.radians(lat)
    lon = np.radians(lon)

    gc = np.arctan(b**2*np.tan(lat)/a**2)
    p = np.sqrt(1./(np.tan(gc)**2/b**2 + 1/a**2))

    z = p * np.tan(gc)
    z += alt * np.sin(lat)
    p += alt * np.cos(lat)
    x = p * np.cos(lon)
    y = p * np.sin(lon)
    return(x, y, z)
