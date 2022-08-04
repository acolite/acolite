## def ged_find
## finds ASTER GED files covering a given limit
## written by Quinten Vanhellemont, RBINS
## 2022-08-03
## modifications: 2022-08-04 (QV) added tile list
##                                switch to usgs opendap server, faster and smaller files
##                                fix latitude in ged_limit

def ged_find(limit, required = False,
             ged_dir = None,
             #ged_url = 'https://e4ftl01.cr.usgs.gov/ASTT/AG100.003/2000.01.01/'
             ged_url = 'https://opendap.cr.usgs.gov/opendap/hyrax/ASTER/ASTT/AG100.003/2000.01.01/'):
    import os, time
    import numpy as np
    import acolite as ac

    if ged_dir is None:
        ged_dir = ac.config['ged_dir'] + '/AG100.003'

    if not os.path.exists(ged_dir): os.makedirs(ged_dir)

    tiles = []
    tilefile = '{}/{}'.format(ac.config['ged_dir'], 'AG100.003_tilelist.txt')
    with open(tilefile, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            tiles.append(line.strip())

    tile_base = 'AG100.v003.{}{}.{}{}.0001.h5'

    ## add one to lat to get correct tile
    ged_limit = [int(np.floor(limit[0]))+1,
                 int(np.floor(limit[1])),
                 int(np.ceil(limit[2]))+1,
                 int(np.ceil(limit[3]))]

    if ged_limit[2] > limit[2]: ged_limit[2]+= -1
    if ged_limit[3] > limit[3]: ged_limit[3]+= -1

    ncol = ged_limit[2] - ged_limit[0] + 1
    nrow = ged_limit[3] - ged_limit[1] + 1

    lat = ged_limit[0]
    lon = ged_limit[1]

    ged_required = []
    for lon in [ged_limit[1] + i for i in range(nrow)]:
        for lat in [ged_limit[0] + j for j in range(ncol)]:

            lat_pf = "-" if lat < 0 else ""
            lon_pf = "-" if lon < 0 else ""

            ged_file = tile_base.format(lat_pf,str(abs(lat)).zfill(2),lon_pf,str(abs(lon)).zfill(3))
            ## test if tile exists
            if ged_file[0:-3] in tiles:
                ged_required.append(ged_file)

    ged_files = []
    for ged_file in ged_required:
        ged_path = '{}/{}'.format(ged_dir, ged_file)
        if 'opendap' in ged_url: ged_path += '.nc4'
        print(ged_path)
        ## try downloading if we don't have the tile
        if not os.path.exists(ged_path):
            ged_local = ac.ged.ged_download(ged_file, ged_url, ged_dir=ged_dir)
            ## add some delay
            #if os.path.exists(ged_local): time.sleep(5)

        if os.path.exists(ged_path): ged_files.append(ged_path)
    return(ged_files)
