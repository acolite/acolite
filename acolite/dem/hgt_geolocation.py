## def hgt_geolocation
## generates geolocation for DEM HGT SRTM files
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-07-17
## last update: 2017-07-18 (QV) added grid keyword
##                2021-04-07 (QV) changed numpy import

def hgt_geolocation(file, grid=True):
    import os
    import numpy as np

    split = os.path.basename(file).split('.')
    lat_s = 1 if split[0][0] == 'N' else -1
    lat_0 = float(split[0][1:3])*lat_s
    lon_s = 1 if split[0][3] == 'E' else -1
    lon_0 = float(split[0][4:7])*lon_s

    dim = (1201,1201)
    step = (1./(dim[0]-1), 1./(dim[1]-1))

    latslice = [lat_0 + i * step[1] for i in range(dim[1])]
    lonslice = [lon_0 + i * step[0] for i in range(dim[0])]

    if grid:
        lon = np.tile(lonslice,dim[1]).reshape(dim)
        lat = np.rot90(np.tile(latslice,dim[0]).reshape(dim))
        return(lon,lat)
    else: return(lonslice,latslice)
