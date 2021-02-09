## def projection_sub
## generic image/projection subset
## written by Quinten Vanhellemont, RBINS
## 2021-02-05
## modifications: 2021-02-09 (QV) added one pixel to sub to correspond to gdalwarp output sizes

def projection_sub(dct, limit, four_corners=True):

    ## input dict has pixel size (x, y)
    ## dimensions (y, x)
    ## scene xrange and yrange
    pixel_size = dct['pixel_size']
    #dimensions = dct['dimensions']
    xdim = dct['xdim']
    ydim = dct['ydim']
    xscene = dct['xrange']
    yscene = dct['yrange']

    if four_corners:
        # West, West, East, East
        # South, North, North, South
        xrange_raw, yrange_raw = dct['p']((limit[1],limit[1],limit[3],limit[3]),
                                          (limit[0],limit[2],limit[2],limit[0]))
        xrange_raw = (min(xrange_raw), max(xrange_raw))
        yrange_raw = (min(yrange_raw), max(yrange_raw))
    else:
        xrange_raw, yrange_raw = dct['p']([limit[1],limit[3]],[limit[0],limit[2]])

    xrange = [xrange_raw[0] - (xrange_raw[0] % pixel_size[0]), xrange_raw[1]+pixel_size[0]-(xrange_raw[1] % pixel_size[0])]
    yrange = [yrange_raw[1]+pixel_size[1]-(yrange_raw[1] % pixel_size[1]), yrange_raw[0] - (yrange_raw[0] % pixel_size[1])]

    if (xrange[1] < xscene[0]) or (xrange[0] > xscene[1]):
        #print('Limits out of scene longitude')
        out_lon = True
    else:
        out_lon = False

    if (yrange[0] < yscene[1]) or (yrange[1] > yscene[0]):
        #print('Limits out of scene latitude')
        out_lat = True
    else:
        out_lat = False

    xrange_region = [x for x in xrange]
    yrange_region = [y for y in yrange]

    xrange = [max((xrange_region[0], xscene[0])), min((xrange_region[1], xscene[1]))]
    yrange = [min((yrange_region[0], yscene[0])), max((yrange_region[1], yscene[1]))]

    xsize = int(xrange[1]-xrange[0])
    ysize = int(yrange[1]-yrange[0])

    xsize_pix = int(xsize/pixel_size[0])
    ysize_pix = int(ysize/pixel_size[1])

    xregion = [int((i - xscene[0])/pixel_size[0]) for i in xrange_region]
    yregion = [int((i - yscene[0])/pixel_size[1]) for i in yrange_region]

    xsize_region = int(xrange_region[1]-xrange_region[0])
    ysize_region = int(yrange_region[1]-yrange_region[0])

    xsize_region_pix = int(xsize_region/pixel_size[0])
    ysize_region_pix = int(ysize_region/pixel_size[1])

    xpos = [int((i - xscene[0])/pixel_size[0]) for i in xrange]
    if xpos[0] < 0: xpos[0]=0
    if xpos[1] >= xdim: xpos[1]=xdim

    ypos = [int((i - yscene[0])/pixel_size[1]) for i in yrange]
    if ypos[0] < 0: ypos[0]=0
    if ypos[1] >= ydim: ypos[1]=ydim

    sub = [xpos[0], ypos[0], xpos[1]-xpos[0]+1, ypos[1]-ypos[0]+1]
    sub_region = [xregion[0], yregion[0], xregion[1]-xregion[0]+1, yregion[1]-yregion[0]+1]

    sub_dct = {'source': dct, 'out_lon': out_lon, 'out_lat': out_lat,
               'sub': sub, 'limit': limit, 'p': dct['p'],
               'xdim': sub[2], 'ydim': sub[3],
               'dimensions': (sub[3], sub[2]),
               'dimensions_xfirst': (sub[2], sub[3]),
               'xrange': xrange, 'xpos': xpos, 'xsize': xsize, 'xsize_pix':xsize_pix,
               'yrange': yrange, 'ypos': ypos, 'ysize': ysize, 'ysize_pix':ysize_pix}

    ## save extended region for merging tiles
    sub_dct['region'] = { 'p': dct['p'], 'sub': sub_region,
                         'xdim': sub_region[2], 'ydim': sub_region[3],
                         'xrange': xrange_region,
                         'xsize_pix':xsize_region_pix,
                         'xregion': xregion, 'xsize':xsize_region,
                         'yrange': yrange_region,
                         'ysize_pix':ysize_region_pix,
                         'yregion': yregion, 'ysize':ysize_region}

    ## copy missing keys from input dict
    for k in dct:
        if k not in sub_dct: sub_dct[k] = dct[k]
        if k not in sub_dct['region']: sub_dct['region'][k] = dct[k]

    return(sub_dct)
