## def geometry
## compute geometry for GOES ABI located at subsatellite point lon_0
##
## based on ac.seviri.vaavza and GOES PUG
## could be merged with ellipsoid and altitude as keywords
##
## Determination of Look Angles to Geostationary Communication Satellites
## https://doi.org/10.1061/(ASCE)0733-9453(1994)120:3(115)
##
## written by Quinten Vanhellemont, RBINS
## 2025-06-05
## modifications:

def geometry(lon, lat, lon_0 = 0,
            r = 42164160, ## GOES altitude,
            a = 6378137, ## GRS80 semi_major_axis
            b = 6356752.31414, ## GRS80 semi_minor_axis
            #a = 6378137, ## WGS84 semi_major_axis
            #b = 6356752.3142, ## WGS84 semi_minor_axis
            ealt = 0,):
    import numpy as np

    ## coordinates in radians
    latr = np.radians(lat)
    lonr = np.radians(lon)
    sinlat = np.sin(latr)
    coslat = np.cos(latr)
    sinlon = np.sin(lonr)
    coslon = np.cos(lonr)

    rlon_0 = np.radians(lon_0)
    coslon0 = np.cos(rlon_0)
    sinlon0 = np.sin(rlon_0)

    # eccentricity
    epsilon = ((a**2 - b**2) / (a**2))**0.5

    # principal radius of curvature in the prime vertical
    N = a / (1 - ((epsilon * sinlat)**2))**0.5

    x_ant = (N + ealt) * coslon * coslat
    y_ant = (N + ealt) * sinlon * coslat
    z_ant = (N * (1 - epsilon**2) + ealt) * sinlat

    x_sat = r * coslon0
    y_sat = r * sinlon0
    z_sat = 0

    x = x_sat - x_ant
    y = y_sat - y_ant
    z = z_sat - z_ant

    e = -1 * sinlon * x + coslon * y
    n = -1 * sinlat * coslon * x - sinlat * sinlon * y + coslat * z
    u = coslat * coslon * x + coslat * sinlon * y + sinlat * z

    ## azimuth and elevation angles
    alpha = np.arctan(e / n)
    nu = np.arctan(u / (e**2 + n**2)**0.5)

    ## southern hemisphere
    alpha[(alpha < 0) & (lat <= 0)] += 2 * np.pi
    ## add pi
    alpha[lat > 0] += np.pi

    ## assume 90 degrees azimuth at sub satellite point
    alpha[(lon == lon_0) & (lat == 0)] = np.pi/2

    vaa = np.degrees(alpha)
    vza = 90 - np.degrees(nu)
    return(vaa, vza)
