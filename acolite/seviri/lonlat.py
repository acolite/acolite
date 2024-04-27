## def lonlat
## compute coordinates for MSG/SEVIRI MSG located at subsatellite point lon_0
## based on msg_nav_utils
##
## written by Quinten Vanhellemont, RBINS
## 2024-04-18
## modifications:

def lonlat(lon_0 = 0.0):
    import numpy as np

    ## earth dimensions and satellite distance
    earth_radius_equator = 6378.169
    earth_radius_pole = 6356.5838
    satellite_distance = 42164.0

    ## image info
    nc, nl = 3712, 3712 ## full disk at non HRV bands
    coff, loff = 1856, 1856 ## full disk centre
    cfac, lfac = 781648343, 781648343 ## factors dependent on pixel resolution in micro rad

    ## projection parameters
    p1 = 1.006803
    p2 = 1737121856

    ## set up pixel dimensions
    column_range = np.arange(0, nc)
    line_range = np.arange(0, nl)
    c, l = np.meshgrid(column_range, line_range)

    ## compute scanning angle
    x = 2**16 * (c - coff) / cfac
    y = 2**16 * (l - loff) / lfac

    cosx = np.cos(x)
    cosy = np.cos(y)
    sinx = np.sin(x)
    siny = np.sin(y)

    sd = np.sqrt((satellite_distance*cosx*cosy)**2 - (cosy**2 + p1*siny**2)*p2)
    sn = (satellite_distance*cosx*cosy - sd) / (cosy**2 + p1*siny**2)
    s1 = satellite_distance - sn*cosx*cosy
    s2 = sn*sinx*cosy
    s3 = -sn*siny
    sxy = np.sqrt(s1**2 + s2**2)

    lon = np.degrees(np.arctan(s2/s1)) + lon_0
    lat = np.degrees(np.arctan(p1*s3/sxy))

    return(lon, lat)
