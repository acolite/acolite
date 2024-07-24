## def lonlat
## compute coordinates for MSG/SEVIRI MSG located at subsatellite point lon_0
## based on msg_nav_utils
##          FCI PUG pdf_mtg_fci_l1_pug.pdf EUM/MTG/USR/13/719113
##
## written by Quinten Vanhellemont, RBINS
## 2024-04-18
## modifications: 2024-07-24 (QV) added FCI, added sensor and ssd (spatial sampling distance) keywords

def lonlat(lon_0 = 0.0, sensor = 'SEVIRI', ssd = 0.5):
    import numpy as np

    if sensor.upper() == 'SEVIRI':
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
    elif sensor.upper() == 'FCI':
        ## earth dimensions and satellite distance (FCI)
        earth_radius_equator = 6378.137
        earth_flattening = 1. / 298.257223563
        geostationary_altitude = 35786.4
        earth_radius_pole = -1*((earth_flattening*earth_radius_equator) - earth_radius_equator)
        satellite_distance = geostationary_altitude + earth_radius_equator
        ## image info
        #nc, nl ## columns and lines
        #coff, loff ## full disk centre
        #cfac, lfac ## factors dependent on pixel resolution in micro rad
        if ssd == 0.5: ## full disk at SSD 0.5 km (VIS 0.6 and NIR2.2)
            nc, nl = 22272, 22272
            coff, loff = 11136.5, 11136.5
            cfac, lfac = [np.round(2**16 / 1.3971788E-05).astype(int)] * 2
        elif ssd == 1.0: ## full disk at SSD 1.0 km (VSWIR and IR 3.8, IR 10.5)
            nc, nl = 11136, 11136
            coff, loff = 5568.5, 5568.5
            cfac, lfac = [np.round(2**16 / 2.7943576E-05).astype(int)] * 2
        elif ssd == 2.0: ## full disk at SSD 2.0 km (TIR)
            nc, nl = 5568, 5568
            coff, loff = 2784.5, 2784.5
            cfac, lfac = [np.round(2**16 / 5.5887153E-05).astype(int)] * 2
        else:
            print('SSD = {} km not configured, use 0.5, 1.0, or 2.0'.format(ssd))
            return
        ## projection parameters (FCI)
        ## p1 == S4 from PUG
        p1 = (earth_radius_equator**2) / (earth_radius_pole**2)
        ## p2 == S5 from PUG
        p2 = np.round((satellite_distance**2) - (earth_radius_equator**2)).astype(int)
    else:
        print('Sensor = {} not configured, use SEVIRI (MSG) or FCI (MTG)')
        return

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
