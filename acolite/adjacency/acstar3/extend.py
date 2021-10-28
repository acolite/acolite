## extend image by mirroring edges
## QV 2021-07-07 from AC R code .extend_scene_v0.9

def extend(d, idp, fill_nan=True):
    import numpy as np
    import acolite as ac
    
    dim = d.shape
    data_ext = np.zeros((dim[0]+idp*2, dim[1]+idp*2))

    ## center
    data_ext[idp:idp+dim[0], idp:idp+dim[1]] = d

    ## top left corner
    data_ext[0:idp, 0:idp] = np.flip(d[0:idp, 0:idp])

    ## left side
    data_ext[idp:idp+dim[0], 0:idp] = np.fliplr(d[:, 0:idp])

    ## bottom left corner
    data_ext[-idp:, 0:idp] = np.flip(d[-idp:, 0:idp])

    ## bottom side
    data_ext[-idp:, idp:-idp] = np.flipud(d[-idp:,:])

    ## bottom right corner
    data_ext[-idp:, -idp:] = np.flip(d[-idp:, -idp:])

    ## right side
    data_ext[idp:idp+dim[0], -idp:] = np.fliplr(d[:, -idp:])

    ## top right corner
    data_ext[0:idp, -idp:] = np.flip(d[0:idp, -idp:])

    ## top edge
    data_ext[0:idp, idp:-idp] = np.flipud(d[0:idp, :])

    ## add mean reflectance for nans
    #if fill_nan: data_ext[np.isnan(data_ext)] = np.nanmean(d)

    ## add closest data for nans
    if fill_nan: data_ext = ac.shared.fillnan(data_ext)

    return(data_ext)
