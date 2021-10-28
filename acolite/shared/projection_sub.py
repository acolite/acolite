## def projection_sub
## generic image/projection subset
## written by Quinten Vanhellemont, RBINS
## 2021-02-05
## modifications: 2021-02-09 (QV) added one pixel to sub to correspond to gdalwarp output sizes
##                2021-10-28 (QV) fix for negative pixel indices and scene grid offset

def projection_sub(dct, limit, four_corners=True, target_pixel_size = None):

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

    if target_pixel_size is None:
        xrange = [xrange_raw[0] - (xrange_raw[0] % pixel_size[0]), xrange_raw[1]+pixel_size[0]-(xrange_raw[1] % pixel_size[0])]
        yrange = [yrange_raw[1]+pixel_size[1]-(yrange_raw[1] % pixel_size[1]), yrange_raw[0] - (yrange_raw[0] % pixel_size[1])]
        ## compute whether we are off from the nominal grid
        x_grid_off = xscene[0]%pixel_size[0], xscene[1]%pixel_size[0]
        y_grid_off = yscene[0]%pixel_size[1], yscene[1]%pixel_size[1]
    else:
        xrange = [xrange_raw[0] - (xrange_raw[0] % target_pixel_size[0]*2), xrange_raw[1]+target_pixel_size[0]*2-(xrange_raw[1] % target_pixel_size[0]*2)]
        yrange = [yrange_raw[1]+target_pixel_size[1]*2-(yrange_raw[1] % target_pixel_size[1]*2), yrange_raw[0] - (yrange_raw[0] % target_pixel_size[1]*2)]
        ## compute whether we are off from the nominal grid
        x_grid_off = xscene[0]%target_pixel_size[0], xscene[1]%target_pixel_size[0]
        y_grid_off = yscene[0]%target_pixel_size[1], yscene[1]%target_pixel_size[1]

    ## remove pixel grid offsets
    xrange[0]+=x_grid_off[0]
    xrange[1]+=x_grid_off[1]
    yrange[0]+=y_grid_off[0]
    yrange[1]+=y_grid_off[1]

    if (xrange[1] < xscene[0]) or (xrange[0] > xscene[1]):
        out_lon = True
    else:
        out_lon = False

    if (yrange[0] < yscene[1]) or (yrange[1] > yscene[0]):
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

    ## if pixel indices are negative remove an extra pixel after converting to int to handle coordinate being "left" or "above" the scene
    xregion = [((i - xscene[0])/pixel_size[0]) for i in xrange_region]
    yregion = [((i - yscene[0])/pixel_size[1]) for i in yrange_region]
    xregion[0] = int(xregion[0]) if xregion[0] > 0 else int(xregion[0])
    xregion[1] = int(xregion[1]) #if xregion[1] > 0 else int(xregion[1])-1
    yregion[0] = int(yregion[0]) if yregion[0] > 0 else int(yregion[0])
    yregion[1] = int(yregion[1]) #if yregion[1] > 0 else int(yregion[1])-1

    xsize_region = int(xrange_region[1]-xrange_region[0])
    ysize_region = int(yrange_region[1]-yrange_region[0])

    xsize_region_pix = int(xsize_region/pixel_size[0])
    ysize_region_pix = int(ysize_region/pixel_size[1])

    #xpos = [int((i - xscene[0])/pixel_size[0]) for i in xrange]
    xpos = [(i - xscene[0])/pixel_size[0] for i in xrange]
    if xpos[0] < 0: xpos[0]=0
    else: xpos[0] = int(xpos[0])
    if xpos[1] >= xdim: xpos[1]=xdim
    else: xpos[1] = int(xpos[1])

    #ypos = [int((i - yscene[0])/pixel_size[1]) for i in yrange]
    ypos = [(i - yscene[0])/pixel_size[1] for i in yrange]
    if ypos[0] < 0: ypos[0]=0
    else: ypos[0] = int(ypos[0])
    if ypos[1] >= ydim: ypos[1]=ydim
    else: ypos[1] = int(ypos[1])

    #sub = [xpos[0], ypos[0], xpos[1]-xpos[0]+1, ypos[1]-ypos[0]+1]
    #sub_region = [xregion[0], yregion[0], xregion[1]-xregion[0]+1, yregion[1]-yregion[0]+1]
    # the adding of a pixel should not be necessary with the above fix for negative coordinates
    sub = [xpos[0], ypos[0], xpos[1]-xpos[0], ypos[1]-ypos[0]]
    sub_region = [xregion[0], yregion[0], xregion[1]-xregion[0], yregion[1]-yregion[0]]

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
