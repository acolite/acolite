## gets satellite elevation and azimuth from x,y,z Earth Centred Earth Fixed
## for a given observation position obs_lat, obs_lon
##
## written by Quinten Vanhellemont, RBINS
## 2025-07-22
## modifications:

def ecef_geometry(sat_x, sat_y, sat_z,
                 obs_lon, obs_lat, obs_alt,
                 a = 6378.137, f = 1/298.257223563, ## WGS84
                 return_elevation = False):
    import numpy as np

    ret = {}
    ## convert lon, lat to radians
    lat = np.radians(obs_lat)
    lon = np.radians(obs_lon)

    ## 1) compute ECEF for obs location
    b = a - (f*a)
    N = a ** 2 / (a ** 2 * np.cos(lat) ** 2 + b ** 2 * np.sin(lat) **2 ) ** 0.5
    x0 = (N + obs_alt) * np.cos(lat) * np.cos(lon)
    y0 = (N + obs_alt) * np.cos(lat) * np.sin(lon)
    z0 = (N * (b / a) ** 2 + obs_alt) * np.sin(lat)

    ## 2) compute East North Up coordinates from difference
    t = np.cos(lon) * (sat_x - x0) + np.sin(lon) * (sat_y - y0)
    E = -np.sin(lon) * (sat_x - x0) + np.cos(lon) * (sat_y - y0)
    U = np.cos(lat) * t + np.sin(lat) * (sat_z - z0)
    N = -np.sin(lat) * t + np.cos(lat) * (sat_z - z0)

    ## 3) compute slant range and zenith, azimuth angles
    r = np.hypot(E, N)
    ret['range'] = np.hypot(r, U)
    zenith = np.pi / 2. - np.arctan2(U, r)
    azimuth = np.arctan2(E, N) % (np.pi * 2.)

    ## convert to degrees
    ret['zenith'] = np.degrees(zenith)
    ret['azimuth'] = np.degrees(azimuth)
    if return_elevation: ret['elevation'] = 90 - ret['zenith']

    return(ret)
