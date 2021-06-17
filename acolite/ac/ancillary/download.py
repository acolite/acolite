## ancillary_download
## gets ancillary data from the ocean data server
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-10-17
## modifications:
##                2018-07-18 (QV) changed acolite import name
##                2018-11-19 (QV) added verbosity option, fixed the download_file script for the GUI
##                2020-03-10 (QV) added autoremove for files that are too small
##                2020-06-23 (QV) changed get url
##                2021-03-01 (QV) simplified for acg renamed from ancillary_download

def download(date, local_dir = None,
                       download=True, override = False, verbosity=0,
                       get_url = "https://oceandata.sci.gsfc.nasa.gov/ob/getfile"):

    import os
    from datetime import datetime, timedelta
    import urllib.request

    import acolite as ac
    if local_dir == None: local_dir=ac.config['met_dir']

    ancillary_files = ac.ac.ancillary.list_files(date)

    local_files = []
    for basefile in ancillary_files:
            yjd = basefile[1:8]
            year = yjd[0:4]
            jday = yjd[4:7]
            url_file = '{}/{}'.format(get_url, basefile)
            local_file = '{}/{}/{}/{}'.format(local_dir,year,jday,basefile)
            local_file_unzipped = local_file.replace('.bz2', '')

            if download:
                ## download file
                ## test if file is large enough
                if os.path.exists(local_file):
                    st = os.stat(local_file)
                    size = st.st_size / (1024 * 1024)
                    if size < 0.05: ## 50Kb
                        if verbosity > 0: print('Removing {} with too small size {:.2f}Mb'.format(os.path.basename(local_file), size))
                        try:
                            os.remove(local_file)
                            if verbosity > 0: print('Deleted {}'.format(local_file))
                        except:
                            if verbosity > 0: print('Could not remove {}'.format(local_file))

                if (os.path.exists(local_file) | os.path.exists(local_file_unzipped)) & (not override):
                    if verbosity > 1: print('File {} exists'.format(basefile))
                    if os.path.exists(local_file_unzipped): local_files.append(local_file_unzipped)
                    elif os.path.exists(local_file): local_files.append(local_file)
                else:
                    if os.path.exists(os.path.dirname(local_file)) is False:
                        os.makedirs(os.path.dirname(local_file))
                    if verbosity > 0: print('Downloading file {}'.format(basefile))
                    try:
                        ac.shared.download_file(url_file, local_file)
                        if verbosity > 0: print('Finished downloading file {}'.format(basefile))
                        local_files.append(local_file)
                    except:
                        print('Downloading file {} failed'.format(basefile))
    return(local_files)
