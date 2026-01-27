## def station_check
## function to convert a station location to a 4 element limit list
## split off from acolite_run
## written by Quinten Vanhellemont, RBINS
## 2026-01-27
## modifications:

def station_check(site_lon, site_lat, box_size, box_units, verbosity = 0):
    import os
    import acolite as ac

    setu = {}
    if verbosity > 0: print('Creating new limit for position {}N, {}E, box size {} {}'.format(site_lat, site_lon, box_size, box_units))
    if box_units[0] in ['k', 'm']:
        if box_units[0] == 'm': box_size /= 1000.
        ## get approximate distance per degree lon/lat
        dlon, dlat = ac.shared.distance_in_ll(site_lat)
        if type(box_size) is list:
            lat_off, lon_off = (box_size[0]/dlat)/2, (box_size[1]/dlon)/2
        else:
           lat_off, lon_off = (box_size/dlat)/2, (box_size/dlon)/2
    elif box_units[0] in ['d']:
        if type(box_size) is list:
            lat_off, lon_off = box_size[0]/2, box_size[1]/2
        else:
            lat_off, lon_off = box_size/2, box_size/2
    else:
        print('station_box_units={} not configured'.format(box_units))
        return
    ## set new limit
    setu['limit'] = [float(site_lat-lat_off), float(site_lon-lon_off), float(site_lat+lat_off), float(site_lon+lon_off)]
    if verbosity > 0: print('New limit: {}'.format(', '.join(['{:}'.format(v) for v in setu['limit']])))
    ## end create limit based on station information
    return(setu)
