## def polygon_check
## function to convert a polygon to a 4 element limit list
## split off from acolite_run
## written by Quinten Vanhellemont, RBINS
## 2026-01-27
## modifications:

def polygon_check(polygon, polygon_limit = True, verbosity = 0):
    import os
    import acolite as ac

    ## empty settings dict
    setu = {}

    if polygon is not None:
        ## is the given polygon a wkt?
        if not os.path.exists(polygon):
            try:
                if (ac.settings['run']['output'] is not None) & ('runid' in ac.settings['run']):
                    polygon_file = '{}/polygon_{}.json'.format(ac.settings['run']['output'], ac.settings['run']['runid'])
                else:
                    polygon_file = '{}/polygon.json'.format(ac.config['scratch_dir'])
                if os.path.exists(polygon_file): os.remove(polygon_file)

                polygon_new = ac.shared.polygon_from_wkt(polygon, file = polygon_file)
                setu['polygon_old'] = '{}'.format(polygon)
                setu['polygon'] = '{}'.format(polygon_new)
                setu['polygon_clip'] = True
            except:
                if verbosity > 1: print('Provided polygon is not a valid WKT polygon')
                if verbosity > 1: print(setu['polygon'])
                setu['polygon'] = None
                pass

        else:
            setu['polygon'] = polygon

        ## read the polygon file
        if os.path.exists(setu['polygon']) & (polygon_limit):
            try:
                limit = ac.shared.polygon_limit(setu['polygon'])
                setu['limit'] = [l for l in limit]
                setu['polygon_clip'] = True
                if verbosity > 1: print('Using limit from polygon envelope: {}'.format(', '.join(['{:}'.format(v) for v in limit])))
            except:
                if verbosity > 1: print('Failed to import polygon {}'.format(setu['polygon']))
        else:
            setu['polygon_clip'] = False

    return(setu)
