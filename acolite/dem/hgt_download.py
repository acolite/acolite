## def hgt_download
## download DEM HGT SRTM files
## written by Quinten Vanhellemont, RBINS
## 2021-04-21
## last update:

def hgt_download(tile,
                 hgt_dir = None, override = False,
                 url_base = 'http://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL3.003/2000.02.11/{}.SRTMGL3.hgt.zip'):

    import os
    import acolite as ac
    if hgt_dir is None: hgt_dir = ac.config['hgt_dir']
    if not os.path.exists(hgt_dir): os.makedirs(hgt_dir)
    
    f_url = url_base.format(tile)
    f_local = '{}/{}'.format(hgt_dir,os.path.basename(f_url))

    if not os.path.exists(f_local):
        try:
            ret = ac.shared.download_file(f_url, f_local)
        except:
            pass

    if os.path.exists(f_local):
        return(f_local)
    else:
        return()
