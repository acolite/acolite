## def copernicus_dem_find
## finds Copernicus DEM TIFs or COGs for given limit
## written by Quinten Vanhellemont, RBINS
## 2022-07-05

def copernicus_dem_find(limit, source = 'copernicus30'):
    import os
    import math
    import acolite as ac

    ## find limits for DEM files
    dem_limit = [int(math.floor(limit[0])),
                 int(math.floor(limit[1])),
                 int(math.ceil(limit[2])),
                 int(math.ceil(limit[3]))]
    if (limit[2] > 0) & (dem_limit[2] > limit[2]): dem_limit[2]+= -1
    if (limit[3] > 0) & (dem_limit[3] > limit[3]): dem_limit[3]+= -1

    ncol = dem_limit[2] - dem_limit[0] + 1
    nrow = dem_limit[3] - dem_limit[1] + 1

    ## find files
    dem_tiles = []
    for lon in [dem_limit[1] + i for i in range(nrow)]:
        for lat in [dem_limit[0] + j for j in range(ncol)]:
            ret = ac.dem.copernicus_dem_download(lon, lat, source)
            if ret is None: continue
            dem_tiles.append(ret)

    return(dem_tiles)
