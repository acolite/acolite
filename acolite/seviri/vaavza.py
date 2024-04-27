## def vaavza
## compute geometry for MSG/SEVIRI MSG located at subsatellite point lon_0
##
## based on  Soler et al. (1994)
## Determination of Look Angles to Geostationary Communication Satellites
## https://doi.org/10.1061/(ASCE)0733-9453(1994)120:3(115)
##
## written by Quinten Vanhellemont, RBINS
## 2024-04-18
## modifications:

def vaavza(lon, lat, lon_0 = 0, ealt = 0):
    import numpy as np

    ## coordinates in radians
    latr = np.radians(lat)
    lonr = np.radians(lon)
    sinlat = np.sin(latr)
    coslat = np.cos(latr)
    sinlon = np.sin(lonr)
    coslon = np.cos(lonr)

    ## SEVIRI altitude
    r = 42164000
    rlon_0 = np.radians(lon_0)
    coslon0 = np.cos(rlon_0)
    sinlon0 = np.sin(rlon_0)

    # Axes of WGS84
    a = 6378137
    b = 6356752.3142

    # eccentricity of WGS84
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
