## def ged_download
## download ASTER GED files
##
## written by Quinten Vanhellemont, RBINS
## 2022-08-03
## modifications: 2022-08-04 (QV) switch to usgs opendap server, faster and smaller files

def ged_download(tile, url_base, ged_dir = None, override = False):
    import zipfile, os
    import acolite as ac

    if ged_dir is None:
        ged_dir = ac.config['ged_dir'] + '/AG100.003'
    if not os.path.exists(ged_dir): os.makedirs(ged_dir)

    for ext in ['.xml', '']:
        f_url = '{}/{}{}'.format(url_base, tile, ext)
        f_local = '{}/{}'.format(ged_dir,os.path.basename(f_url))
        if 'opendap.cr.usgs.gov' in f_url:
            if ext == '.xml': continue
            f_url+='.nc4?'
            f_local+='.nc4'

        print(f_url)
        print(f_local)

        if not os.path.exists(f_local):
            try:
                ret = ac.shared.download_file(f_url, f_local)
                if os.path.exists(f_local):
                    print('Downloaded {}'.format(f_local))
            except BaseException as err:
                print("Download error {}, {}".format(err, type(err)))
                print('Downloading {} failed'.format(f_local))
                pass
        else:
            print('{} exists'.format(f_local))

    if os.path.exists(f_local):
        return(f_local)
    else:
        return()
