## def eci_geometry
## gets satellite elevation and azimuth from x,y,z Earth Centred Inertial coordinates
## for a given observation position obs_lat, obs_lon
##
## written by Quinten Vanhellemont, RBINS
## based on notes from Satellite Times by Dr. T.S. Kelso (https://celestrak.com/columns/)
## 2022-02-09
## modifications: 2022-02-10 (QV) fixed coordinate system rotation
##                2022-02-12 (QV) support obs arrays
##                2025-05-20 (QV) split spherical/oblate options, added variables as keywords
##                2025-05-21 (QV) added SSP
##                2025-05-22 (QV) split off GMST and SSP

def eci_geometry(sat_x, sat_y, sat_z,
                 obs_lon, obs_lat, obs_alt,
                 time,
                 omega_e = 7.29211510e10-5, ## Earth rotation speed (radians/second)
                 Re = 6378.135, ## Earth radius (km)
                 #a = 6378.135, f = 1/298.26, ## WGS72
                 a = 6378.137, f = 1/298.257223563, ## WGS84
                 return_elevation = False, spherical = False):
    import numpy as np
    import datetime, dateutil.parser
    import acolite as ac

    ## get theta_g
    theta_g = ac.shared.position.gmst(time)

    ## convert to arrays
    obs_lat = np.atleast_1d(obs_lat)
    obs_lon = np.atleast_1d(obs_lon)
    obs_alt = np.atleast_1d(obs_alt)

    ## convert to radians
    lat = np.radians(obs_lat)
    lon = np.radians(obs_lon)

    ## compute Greenwich Mean Sidereal Time
    theta = theta_g + lon
    theta = np.mod(theta, np.pi * 2)

    ## dict to store results
    ret = {}

    ## compute observer ECI coordinates
    ## for spherical earth
    if spherical:
        obs_z = (Re + obs_alt) * np.sin(lat)
        R = Re * np.cos(lat)
        obs_x = R * np.cos(theta)
        obs_y = R * np.sin(theta)

        ## range vector
        rx, ry, rz = sat_x - obs_x, sat_y - obs_y, sat_z - obs_z

        ## rotate coordinate system to topocentric-horizon
        rS = np.sin(lat) * np.cos(theta) * rx + np.sin(lat) * np.sin(theta) * ry - np.cos(lat) * rz
        rE = -np.sin(theta) * rx + np.cos(theta)*ry
        rZ = np.cos(lat) * np.cos(theta) * rx + np.cos(lat) *  np.sin(theta) * ry + np.sin(lat) * rz

        ## range
        r = (rS**2 + rE**2 + rZ**2)**0.5

        ## Elevation
        El = np.arcsin(rZ / r)

        ## Azimuth
        Az = np.arctan(-rE / rS)
        Az[rS>0] += np.pi
        Az[Az<0] += np.pi * 2

        del rx, ry, rz, r
        del rS, rE, rZ

        ret['type'] = 'spherical'
        ret['zenith'] = 90-np.degrees(El)
        ret['azimuth'] =  np.degrees(Az)
        if return_elevation: ret['elevation'] = np.degrees(El)
        del Az, El

    ## for oblate earth
    else:
        b = a * (1 - f)
        C = 1 / (1+f*(f-2)*(np.sin(lat)**2))**0.5
        S = ((1-f)**2)*C

        obs_xprime = a * C * np.cos(lat) * np.cos(theta)
        obs_yprime = a * C * np.cos(lat) * np.sin(theta)
        obs_zprime = a * S * np.sin(lat)

        ## range vector
        rxprime, ryprime, rzprime = sat_x - obs_xprime, sat_y - obs_yprime, sat_z - obs_zprime

        ## rotate coordinate system to topocentric-horizon
        rSprime = np.sin(lat)*np.cos(theta)*rxprime + np.sin(lat)*np.sin(theta)*ryprime - np.cos(lat)*rzprime
        rEprime = -np.sin(theta) * rxprime + np.cos(theta)*ryprime
        rZprime = np.cos(lat)*np.cos(theta)*rxprime + np.cos(lat)*np.sin(theta)*ryprime + np.sin(lat)*rzprime

        ## range
        rprime = (rSprime**2 + rEprime**2 + rZprime**2)**0.5

        ## Elevation
        Elprime = np.arcsin(rZprime / rprime)

        ## Azimuth
        Azprime = np.arctan(-rEprime / rSprime)
        Azprime[rSprime>0] += np.pi
        Azprime[Azprime<0] += np.pi * 2

        del rxprime, ryprime, rzprime, rprime
        del rSprime, rEprime, rZprime

        ret['type'] = 'oblate'
        ret['zenith'] = 90-np.degrees(Elprime)
        ret['azimuth'] =  np.degrees(Azprime)
        if return_elevation: ret['elevation'] = np.degrees(Elprime)
        del Azprime, Elprime

    ## return single values if 1D
    for k in ret:
        if k == 'type': continue
        if ret[k].shape == (1,): ret[k] = ret[k][0]

    return(ret)
