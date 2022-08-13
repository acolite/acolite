## def image_resolution
## estimates resolution for unprojected image
## written by Quinten Vanhellemont, RBINS
## 2022-08-14
## modifications:

def image_resolution(input, ellipsoid = 'WGS84', decimals = 0):
    import acolite as ac
    import pyproj
    import numpy as np

    ## get corner coordinates
    if type(input) is tuple:
        limit, (xdim, ydim), corners = input
    else:
        limit, (xdim, ydim), corners = ac.shared.image_extent(input)

    # estimate resolution from great circle distance between corner coordinates
    g = pyproj.Geod(ellps=ellipsoid)

    ## estimate resolution for each side of the image
    upper_res = g.inv(corners['UL']['lon'], corners['UL']['lat'], corners['UR']['lon'], corners['UR']['lat'])[2]/xdim
    lower_res = g.inv(corners['LL']['lon'], corners['LL']['lat'], corners['LR']['lon'], corners['LR']['lat'])[2]/xdim
    right_res = g.inv(corners['UR']['lon'], corners['UR']['lat'], corners['LR']['lon'], corners['LR']['lat'])[2]/ydim
    left_res = g.inv(corners['UL']['lon'], corners['UL']['lat'], corners['LL']['lon'], corners['LL']['lat'])[2]/ydim

    ## compute average from sides of image
    mean_res = np.mean((upper_res,lower_res,right_res,left_res))

    ## round to # decimals
    resolution = np.round(mean_res, decimals=decimals)

    return(resolution)
