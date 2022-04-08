## def download_file
## download_file with authorisation option
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-14
## modifications: 2018-11-19 (QV) added verbosity option, removed the parallel download
##                2020-04-24 (QV) earthdata login
##                2021-05-31 (QV) added local directory creation
##                2022-04-08 (QV) download to scratch directory first

def download_file(url, file, auth = None, session = None,
                    parallel = False, verbosity = 0, verify_ssl = True):

    import requests, time, os, shutil
    import acolite as ac

    file_path = os.path.abspath(file)
    file_dir = os.path.dirname(file_path)
    if not os.path.exists(file_dir): os.makedirs(file_dir)

    ## first download to temp location
    bn = os.path.basename(file_path)
    temp_file = '{}/{}'.format(ac.config['scratch_dir'], bn)
    if os.path.exists(temp_file): os.remove(temp_file)

    start = time.time()

    if ('https://oceandata.sci.gsfc.nasa.gov/' in url) or ('MEASURES/SRTMGL3.003/2000.02.11' in url):
        if ('EARTHDATA_u' in os.environ) & ('EARTHDATA_p' in os.environ):
            username = os.environ['EARTHDATA_u']
            password = os.environ['EARTHDATA_p']
            auth = (username, password)
        else:
            print('EARTHDATA user name and password required for download of {}'.format(url))
            return()

    with requests.Session() as session:
            r1 = session.request('get', url, verify=verify_ssl)
            r = session.get(r1.url, auth=auth, verify=verify_ssl)

            if (r.ok):
                with open(temp_file, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1024*1024):
                        if chunk: # filter out keep-alive new chunks
                            f.write(chunk)
            else:
                if verbosity > 2: print(r.text)
                raise Exception("File download failed")

    ## copy temp file
    if os.path.exists(temp_file):
        shutil.copyfile(temp_file, file_path)
        os.remove(temp_file)

    if verbosity > 1:
        print("Downloaded {}, elapsed Time: {:.1f}s".format(url, time.time() - start))
