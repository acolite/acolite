## def hgt_find
## finds DEM HGT SRTM files covering a given limit
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-07-17
##                2019-04-24 (QV) added support for zip files
##                2021-04-07 (QV) changed numpy import
##                2021-04-21 (QV) added download

def hgt_find(limit, required=False, hgt_dir=None, hgt_ext='.SRTMGL3.hgt.gz'):
    import os
    import numpy as np
    import acolite as ac

    hgt_limit = [int(np.floor(limit[0])),
                 int(np.floor(limit[1])),
                 int(np.ceil(limit[2])),
                 int(np.ceil(limit[3]))]

    if hgt_limit[2] > limit[2]: hgt_limit[2]+= -1
    if hgt_limit[3] > limit[3]: hgt_limit[3]+= -1

    ncol = hgt_limit[2] - hgt_limit[0] + 1
    nrow = hgt_limit[3] - hgt_limit[1] + 1

    lat = hgt_limit[0]
    lon = hgt_limit[1]

    hgt_required = []

    for lon in [hgt_limit[1] + i for i in range(nrow)]:
        for lat in [hgt_limit[0] + j for j in range(ncol)]:

            lat_pf = "S" if lat < 0 else "N"
            lon_pf = "W" if lon < 0 else "E"

            hgt_file = '{}{}{}{}'.format(lat_pf,str(abs(lat)).zfill(2),lon_pf,str(abs(lon)).zfill(3))
            hgt_required.append(hgt_file)

    hgt_files = []

    for hgt_file in hgt_required:
        hgt_path = '{}/{}{}'.format(hgt_dir, hgt_file, hgt_ext)
        ## try downloading if we don't have the tile
        ## note that some tiles do not exist
        if not os.path.exists(hgt_path): ac.dem.hgt_download(hgt_file, hgt_dir=hgt_dir)

        if os.path.exists(hgt_path): hgt_files.append(hgt_path)
        else:
            hgt_path = '{}/{}{}'.format(hgt_dir, hgt_file, hgt_ext.replace('.gz','.zip'))
            if os.path.exists(hgt_path): hgt_files.append(hgt_path)
            else:
                hgt_path = '{}/{}{}'.format(hgt_dir, hgt_file, hgt_ext.replace('.SRTMGL3.hgt.gz','.hgt'))
                if os.path.exists(hgt_path): hgt_files.append(hgt_path)

    if required:
        return(hgt_files, hgt_required)
    else:
        return(hgt_files)
