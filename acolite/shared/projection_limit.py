## def projection_limit
## generic image/projection subset to limit
## written by Quinten Vanhellemont, RBINS
## 2021-02-23
## modifications:

def projection_limit(dct, sub, four_corners=False, add_half_pixel=False):

    ## input dict has pixel size (x, y)
    ## dimensions (y, x)
    ## scene xrange and yrange
    pixel_size = dct['pixel_size']
    xdim = dct['xdim']
    ydim = dct['ydim']
    xscene = dct['xrange']
    yscene = dct['yrange']

    ## sub[0] is xmin in image, sub[2] is xsize in image
    ## sub[1] is ymin in image, sub[3] is ysize in image
    xlim0 = xscene[0] + (sub[0]) * pixel_size[0]
    xlim1 = xscene[0] + (sub[0]+sub[2]-2) * pixel_size[0] ## added -2 to match projection_sub output
    ylim0 = yscene[0] + (sub[1]-1) * pixel_size[1]        ## added -1 to match projection_sub output
    ylim1 = yscene[0] + (sub[1]+sub[3]-1) * pixel_size[1] ## added -1 to match projection_sub output

    if add_half_pixel:
        xlim0 += dct['pixel_size'][0]/2
        xlim1 += dct['pixel_size'][0]/2
        ylim0 += dct['pixel_size'][1]/2
        ylim1 += dct['pixel_size'][1]/2

    if four_corners:
        lon_r, lat_r = dct['p']((xlim0,xlim0,xlim1,xlim1),
                                (ylim0,ylim1,ylim0,ylim1), inverse=True)
    else:
        lon_r, lat_r = dct['p']((xlim0,xlim1),
                                (ylim0,ylim1), inverse=True)

    limit = [min(lat_r), min(lon_r), max(lat_r), max(lon_r)]
    return(limit)
