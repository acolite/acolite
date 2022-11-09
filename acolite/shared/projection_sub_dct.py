## def projection_sub_dct
## generic image/projection subset when sub is given, based on projection_dct
## written by Quinten Vanhellemont, RBINS
## 2022-11-09

def projection_sub_dct(dct, sub):
    import numpy as np

    ## input dict has pixel size (x, y)
    ## dimensions (y, x)
    ## scene xrange and yrange
    pixel_size = dct['pixel_size']
    xdim = dct['xdim']
    ydim = dct['ydim']
    xscene = dct['xrange']
    yscene = dct['yrange']

    ## compute x and y range for given sub
    xrange_new = [dct['xrange'][0] + (sub[0] * dct['pixel_size'][0]),
                  dct['xrange'][0] + ((sub[0]+sub[2]) * dct['pixel_size'][0])]
    yrange_new = [dct['yrange'][0] + (sub[1] * dct['pixel_size'][1]),
                  dct['yrange'][0] + ((sub[1]+sub[3]) * dct['pixel_size'][1])]

    ## run through the same processing as projection_sub
    ## even though this should not be necessary
    xrange_region = [x for x in xrange_new]
    yrange_region = [y for y in yrange_new]

    xrange = [max((xrange_region[0], xscene[0])), min((xrange_region[1], xscene[1]))]
    yrange = [min((yrange_region[0], yscene[0])), max((yrange_region[1], yscene[1]))]

    xsize = int(np.round(xrange[1]-xrange[0]))
    ysize = int(np.round(yrange[1]-yrange[0]))

    xsize_pix = int(np.round(xsize/pixel_size[0]))
    ysize_pix = int(np.round(ysize/pixel_size[1]))

    ## if pixel indices are negative remove an extra pixel after converting to int to handle coordinate being "left" or "above" the scene
    xregion = [((i - xscene[0])/pixel_size[0]) for i in xrange_region]
    yregion = [((i - yscene[0])/pixel_size[1]) for i in yrange_region]
    xregion[0] = int(xregion[0]) if xregion[0] > 0 else int(xregion[0])
    xregion[1] = int(xregion[1])
    yregion[0] = int(yregion[0]) if yregion[0] > 0 else int(yregion[0])
    yregion[1] = int(yregion[1])

    xsize_region = int(np.round(xrange_region[1]-xrange_region[0]))
    ysize_region = int(np.round(yrange_region[1]-yrange_region[0]))

    xsize_region_pix = int(np.round(xsize_region/pixel_size[0]))
    ysize_region_pix = int(np.round(ysize_region/pixel_size[1]))

    xpos = [np.round((i - xscene[0])/pixel_size[0]) for i in xrange]
    if xpos[0] < 0: xpos[0]=0
    else: xpos[0] = int(xpos[0])
    if xpos[1] >= xdim: xpos[1]=xdim
    else: xpos[1] = int(xpos[1])

    ypos = [np.round((i - yscene[0])/pixel_size[1])for i in yrange]
    if ypos[0] < 0: ypos[0]=0
    else: ypos[0] = int(ypos[0])
    if ypos[1] >= ydim: ypos[1]=ydim
    else: ypos[1] = int(ypos[1])

    ## make new sub and sub_dct
    sub = [xpos[0], ypos[0], xpos[1]-xpos[0], ypos[1]-ypos[0]]
    sub_region = [xregion[0], yregion[0], xregion[1]-xregion[0], yregion[1]-yregion[0]]
    sub_dct = {'source': dct, #'out_lon': out_lon, 'out_lat': out_lat,
               'sub': sub, #'limit': limit, 
               'p': dct['p'],
               'xdim': sub[2], 'ydim': sub[3],
               'dimensions': (sub[3], sub[2]),
               'dimensions_xfirst': (sub[2], sub[3]),
               'xrange': xrange, 'xpos': xpos, 'xsize': xsize, 'xsize_pix':xsize_pix,
               'yrange': yrange, 'ypos': ypos, 'ysize': ysize, 'ysize_pix':ysize_pix}

    ## copy missing keys from input dict
    for k in dct:
        if k not in sub_dct: sub_dct[k] = dct[k]
    return(sub_dct)
