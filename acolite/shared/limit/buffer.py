## def limit.buffer
## adds a buffer to a 4 element limit list
## split off from acolite_run
## written by Quinten Vanhellemont, RBINS
## 2026-01-27
## modifications:

def buffer(limit, limit_buffer, limit_buffer_units, verbosity = 0):
    import os
    import acolite as ac

    ## updated settings
    setu = {}

    ## check limit and add buffer if needed
    if limit is not None:
        if len(limit) != 4:
            print('ROI limit should be four elements in decimal degrees: limit=S,W,N,E')
            print('Provided in the settings:', limit)
            return

        ## add limit buffer
        if (limit_buffer is not None) & (limit_buffer_units is not None):
            if verbosity > 1: print('Applying limit buffer of {} {}'.format(limit_buffer, limit_buffer_units))
            setu['limit_old'] = [l for l in limit]
            setu['limit_buffer'] = limit_buffer

            if limit_buffer_units[0].lower() == 'd':
                limit_factor = 1.0, 1.0
            elif limit_buffer_units[0].lower() in ['m','k']:
                mean_lat = (limit[0] + limit[2]) / 2.
                dlon, dlat = ac.shared.distance_in_ll(lat = mean_lat)
                limit_factor = 1/dlat, 1/dlon
                if limit_buffer_units[0].lower() == 'm':
                    limit_factor = limit_factor[0]/1000, limit_factor[1]/1000
            else:
                print('limit_buffer_units={} not configured'.format(limit_buffer_units))
                return

            setu['limit_factor'] = float(limit_factor[0]), float(limit_factor[1])

            ## compute limit buffer
            setu['limit_buffer_degrees'] = setu['limit_buffer'] * setu['limit_factor'][0], \
                                           setu['limit_buffer'] * setu['limit_factor'][1]
            setu['limit'] = [setu['limit_old'][0] - setu['limit_buffer_degrees'][0], \
                             setu['limit_old'][1] - setu['limit_buffer_degrees'][1], \
                             setu['limit_old'][2] + setu['limit_buffer_degrees'][0], \
                             setu['limit_old'][3] + setu['limit_buffer_degrees'][1]]

            if verbosity > 1:
                print('Old limit: {}'.format(', '.join(['{:}'.format(v) for v in setu['limit_old']])))
                print('New limit: {}'.format(', '.join(['{:}'.format(v) for v in setu['limit']])))
    return(setu)
