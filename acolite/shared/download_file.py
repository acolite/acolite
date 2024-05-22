## def download_file
## download_file with authorisation option
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-14
## modifications: 2018-11-19 (QV) added verbosity option, removed the parallel download
##                2020-04-24 (QV) earthdata login
##                2021-05-31 (QV) added local directory creation
##                2022-04-08 (QV) download to scratch directory first
##                2022-07-07 (QV) added SRTM1 DEM
##                2022-08-04 (QV) added GED and retry option
##                2022-08-17 (QV) added .netrc auth, simplified url checks for earthdata
##                2024-05-01 (QV) added earthdatacloud.nasa.gov check for earthdata
##                2024-05-22 (QV) use EARTHDATA_urls from config

def download_file(url, file, auth = None, session = None,
                    parallel = False, verbosity = 0, verify_ssl = True, retry = 1):

    import requests, time, os, shutil, netrc
    import acolite as ac

    file_path = os.path.abspath(file)
    file_dir = os.path.dirname(file_path)
    if not os.path.exists(file_dir): os.makedirs(file_dir)

    ## first download to temp location
    bn = os.path.basename(file_path)
    temp_file = '{}/{}'.format(ac.config['scratch_dir'], bn)
    if os.path.exists(temp_file): os.remove(temp_file)

    start = time.time()

    if any([u in url for u in ac.config['EARTHDATA_urls']]):
        ## try to get auth from netrc
        if (auth is None):
            try:
                nr = netrc.netrc()
                ret = nr.authenticators('earthdata')
                if ret is not None:
                    login, account, password = ret
                    login = login.strip('"')
                    password = password.strip('"')
                    auth = (login, password)
            except:
                pass

        if (auth is None) & ('EARTHDATA_u' in os.environ) & ('EARTHDATA_p' in os.environ):
            username = os.environ['EARTHDATA_u']
            password = os.environ['EARTHDATA_p']
            auth = (username, password)
        if (auth is None):
            print('EARTHDATA user name and password required for download of {}'.format(url))
            return()

    with requests.Session() as session:
            r1 = session.request('get', url, verify=verify_ssl)
            r = session.get(r1.url, auth=auth, verify=verify_ssl)
            time.sleep(1)

            if (r.ok):
                with open(temp_file, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1024*1024):
                        if chunk: # filter out keep-alive new chunks
                            f.write(chunk)
            else:
                if verbosity > 2: print(r.text)
                if retry > 0:
                    retry -= 1
                    print('Retrying...')
                    ac.shared.download_file(url, file, auth = auth, session = session,
                                                parallel = parallel, verbosity = verbosity, verify_ssl = verify_ssl, retry = retry)
                if not os.path.exists(file):
                    raise Exception("File download failed {}".format(r.text))

    ## copy temp file
    if os.path.exists(temp_file):
        shutil.copyfile(temp_file, file_path)
        os.remove(temp_file)

    if verbosity > 1:
        print("Downloaded {}, elapsed Time: {:.1f}s".format(url, time.time() - start))
