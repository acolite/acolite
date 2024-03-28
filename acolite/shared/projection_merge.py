## def projection_merge
## merge list of dct projections
##
## written by Quinten Vanhellemont, RBINS
## 2024-03-28
## modifications: 2024-03-28 (QV) split from multi_tile_extent

def projection_merge(dct_list, dct = None):
    import numpy as np

    for dct_ in dct_list:
        if dct is None:
            dct = {k: dct_[k] for k in dct_}
        else:
            if dct_['epsg'] == dct['epsg']:
                dct['xrange'][0] = np.min((dct['xrange'][0], dct_['xrange'][0]))
                dct['xrange'][1] = np.max((dct['xrange'][1], dct_['xrange'][1]))
                dct['yrange'][1] = np.min((dct['yrange'][1], dct_['yrange'][1]))
                dct['yrange'][0] = np.max((dct['yrange'][0], dct_['yrange'][0]))
            else:
                ### if the prj does not match, project current scene bounds to lat/lon
                lonr, latr = dct_['p']([dct_['xrange'][0], dct_['xrange'][0],
                                        dct_['xrange'][1], dct_['xrange'][1]],
                                       [dct_['yrange'][1], dct_['yrange'][0],
                                        dct_['yrange'][1], dct_['yrange'][0]], inverse=True)

                ### then to target projection
                xrange_raw, yrange_raw = dct['p'](lonr, latr)
                xrange_raw = [np.min(xrange_raw), np.max(xrange_raw)]
                yrange_raw = [np.min(yrange_raw), np.max(yrange_raw)]

                ## fix to nearest full pixel
                xrange = [xrange_raw[0] - (xrange_raw[0] % dct['pixel_size'][0]), xrange_raw[1]+dct['pixel_size'][0]-(xrange_raw[1] % dct['pixel_size'][0])]
                yrange = [yrange_raw[1]+dct['pixel_size'][1]-(yrange_raw[1] % dct['pixel_size'][1]), yrange_raw[0] - (yrange_raw[0] % dct['pixel_size'][1])]

                ## update dct
                dct['xrange'][0] = np.min((dct['xrange'][0], xrange[0]))
                dct['xrange'][1] = np.max((dct['xrange'][1], xrange[1]))
                dct['yrange'][1] = np.min((dct['yrange'][1], yrange[1]))
                dct['yrange'][0] = np.max((dct['yrange'][0], yrange[0]))

            ## update dimensions
            dct['xdim'] = int((dct['xrange'][1]-dct['xrange'][0]) / dct['pixel_size'][0])
            dct['ydim'] = int((dct['yrange'][1]-dct['yrange'][0]) / dct['pixel_size'][1])
            dct['dimensions'] = dct['xdim'], dct['ydim']
    return(dct)
