## def hgt_lonlat
## gets SRTM DEM for given lon/lat arrays
## written by Quinten Vanhellemont, RBINS
## QV 2017-07-17
## modifications: 2018-12-04 QV added printing of url to download (USGS server requires login so presently not automated)
##                2021-04-07 (QV) added to generic acolite
##                2021-04-21 (QV) removed return if tiles are missing (this is also possible since hgt_find does not know which tiles exist)
##                2022-07-07 (QV) added SRTM1 DEM

def hgt_lonlat(lon1, lat1, nearest=True, hgt_dir=None, source = 'srtmgl3'):

    import os
    import acolite as ac
    from scipy import interpolate
    import numpy as np

    if hgt_dir is None: hgt_dir = ac.config['hgt_dir']

    if source == 'srtmgl3':
        url_base = 'https://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL3.003/2000.02.11/{}.SRTMGL3.hgt.zip'
        hgt_dir+='/SRTMGL3'
    elif source == 'srtmgl1':
        url_base = 'https://e4ftl01.cr.usgs.gov/DP133/SRTM/SRTMGL1.003/2000.02.11/{}.SRTMGL1.hgt.zip'
        hgt_dir+='/SRTMGL1'

    ## find dem files
    limit=[0,0,0,0]
    if (type(lon1) is float) | (type(lon1) is int):
        limit[1]=lon1
        limit[3]=lon1
    else:
        limit[1]=lon1.min()
        limit[3]=lon1.max()

    if (type(lat1) == float) | (type(lat1) == int):
        limit[0]=lat1
        limit[2]=lat1
    else:
        limit[0]=lat1.min()
        limit[2]=lat1.max()

    hgt_files, hgt_required = ac.dem.hgt_find(limit, required=True, hgt_dir=hgt_dir, hgt_url=url_base)

    #if len(hgt_files) != len(hgt_required):
    #    print('DEM files not found in {}'.format(hgt_dir))
    #    for f in hgt_required:
    #        if f in hgt_files: continue
    #        f_url = url_base.format(f)
    #        f_local = '{}/{}'.format(hgt_dir,os.path.basename(f_url))
    #        print('Please download {} to {}'.format(f_url,f_local))
    #    #return(0)

    dem = np.asarray(0.0)

    ## run through dem files and reproject data to target lat,lon
    for i, hgt_file in enumerate(hgt_files):
        ## read hgt data and geolocation
        hgt = ac.dem.hgt_read(hgt_file)

        if (type(lon1) is float) & (type(lat1) is float):
            lon0,lat0 = ac.dem.hgt_geolocation(hgt_file, grid=False)
            hgtip = interpolate.RectBivariateSpline(lon0,lat0, hgt)
            result = hgtip(lon1,lat1)
        else:
            lon0,lat0 = ac.dem.hgt_geolocation(hgt_file, grid=True)
            ## reproject
            result = ac.shared.reproject2(hgt, lon0, lat0, lon1, lat1, nearest=nearest)

        ## make output
        if i == 0:
            dem = result
        else:
            ## copy new results where != 0 (pyresample fill does not seem to be working)
            dem[result != 0] = result[result != 0]

    return(dem)
