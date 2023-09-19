## def download
## download urls combinations from CDSE
## written by Quinten Vanhellemont, RBINS
## 2023-09-12
## modifications:

def download(urls, output = None, auth = None, auth_url = None,
                  extract_zip = True, remove_zip = True, override = False, verbosity = 1):
    import os, requests, netrc
    import acolite as ac

    ## get authentication URL from config
    if auth_url is None: auth_url = ac.config['CDSE_auth']

    ## get auth for download
    if auth is None:

        try:
            ## get auth from netrc file
            nr = netrc.netrc()
            ret = nr.authenticators('cdse')
            if ret is not None:
                login, account, password = ret
                login = login.strip('"')
                password = password.strip('"')
                auth = (login, password)
        except:
            pass

        ## get auth from environment
        if auth is None:
            if ('CDSE_u' in os.environ) & \
               ('CDSE_p' in os.environ):
                auth = (os.environ['CDSE_u'], os.environ['CDSE_p'])

        ## get auth from config
        if auth is None:
            if (ac.config['CDSE_u'] != '') & \
               (ac.config['CDSE_p'] != ''):
                auth = (ac.config['CDSE_u'], ac.config['CDSE_p'])

    if auth is None:
        print('Could not determine CDSE credentials CDSE_u and CDSE_p.')
        print('Please add them to your .netrc file, environment varialbles, or ACOLITE config file.')
        return()

    ## get access token
    data = {"client_id": "cdse-public", "grant_type":
            "password","username": auth[0], "password": auth[1]}
    response = requests.post(auth_url, data=data, verify=True, allow_redirects=False)
    access_token = response.json()['access_token']
    if verbosity > 1: print(access_token)

    if output is None:
        cwd = os.getcwd()
        if verbosity > 0: print('No output directory given, will download to current working directory: {}'.format(cwd))
        output = '{}'.format(cwd)

    ## create path
    if not os.path.exists(output): os.makedirs(output)

    ## set up download session
    session = requests.Session()
    session.headers["Authorization"] = f"Bearer {access_token}"

    ## still convert to list if only one url or scene given
    if type(urls) is str: urls = [urls]

    ## run through urls and download the file
    if verbosity > 0: print('Downloading {} scenes'.format(len(urls)))
    lfiles, zfiles = [], []
    for ui, url in enumerate(urls):
        ## find out scene name
        response = requests.get(url.strip('/$value'))
        scene_atts = response.json()
        scene = scene_atts['value'][0]['Name']

        ## local files
        lfile = '{}/{}'.format(output, scene)
        zfile = '{}/{}.zip'.format(output, scene[0:scene.find('.')])

        ## download if we don't have the scene
        if (override | (not os.path.exists(lfile)) & (not os.path.exists(zfile))):
            ## try url
            print('Downloading {}'.format(scene))
            response = session.get(url, allow_redirects=False)
            ## follow redirects
            while response.status_code in (301, 302, 303, 307):
                url = response.headers['Location']
                response = session.get(url, allow_redirects=False)

            ## download file
            dl = session.get(url, verify=False, allow_redirects=True)
            print('Writing file to {}'.format(zfile))
            if (dl.ok):
                with open(zfile, 'wb') as p:
                    for chunk in dl.iter_content(chunk_size=1024*1024):
                        if chunk: # filter out keep-alive new chunks
                            p.write(chunk)
            else:
                print('An error occurred trying to download.')

        else:
            print('Local copy of {} exists'.format(scene))

        ## extract zip file
        if extract_zip:
            if os.path.exists(zfile):
                print('Extracting {}'.format(zfile))
                ac.shared.extract_bundle(zfile, targ_bundle=lfile)
                print('Wrote {}'.format(lfile))

            ## remove downloaded zip file after extraction
            if remove_zip:
                if os.path.exists(zfile):
                    print('Deleting {}'.format(zfile))
                    os.remove(zfile)
                    print('Deleted {}'.format(zfile))

        ## list local paths
        if extract_zip:
            if os.path.exists(lfile): lfiles.append(lfile)
        else:
            if os.path.exists(zfile): lfiles.append(zfile)
    return(lfiles)
