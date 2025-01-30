## def projection_warp_to
## make warp to tuple from projection dict
## written by Quinten Vanhellemont, RBINS
## 2025-01-25
## modifications:

def projection_warp_to(inp, res_method = 'bilinear'):
    import os
    import acolite as ac

    if type(inp) is str:
        if os.path.exists(inp): dct = ac.shared.projection_read(inp)
    elif type(inp) is dict:
        dct = {k: inp[k] for k in inp}
    else:
        dct = None

    if dct is None:
        print('Could not set up warp_to projection from dct')
        print(dct)
        return

    projection = None
    if 'projection' in dct:
        projection = '{}'.format(dct['projection'])
    elif 'epsg' in dct:
        projection = '{}'.format(dct['epsg'])
    else:
        print('Could not determine projection from dct')
        print(dct)
        return

    xyr = [min(dct['xrange']),
           min(dct['yrange']),
           max(dct['xrange']),
           max(dct['yrange']),
           projection]

    warp_to = (projection, xyr,
                dct['pixel_size'][0],
                dct['pixel_size'][1], res_method)

    return(warp_to)
