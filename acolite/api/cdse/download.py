## def download
## download urls combinations from CDSE
## written by Quinten Vanhellemont, RBINS
## 2023-09-12
## modifications: 2023-09-19 (QV) retrieve access token per url to avoid time outs
##                2023-10-22 (QV) added optional scenes list
##                2024-04-27 (QV) moved to acolite.api

def download(urls, scenes = [], output = None, auth = None, auth_url = None, netrc_machine = 'cdse',
                  extract_zip = True, remove_zip = True, override = False, verbosity = 1):
    import os, requests, netrc, time
    import acolite as ac

    ## get authentication URL from config
    if auth_url is None: auth_url = ac.config['CDSE_auth']

    ## get credentials
    if auth is None: auth = ac.shared.auth(netrc_machine)

    if auth is None:
        print('Could not determine CDSE credentials CDSE_u and CDSE_p.')
        print('Please add them to your .netrc file, environment variables, or ACOLITE config file.')
        return()

    if output is None:
        cwd = os.getcwd()
        if verbosity > 0: print('No output directory given, will download to current working directory: {}'.format(cwd))
        output = '{}'.format(cwd)

    ## create path
    if not os.path.exists(output): os.makedirs(output)

    ## still convert to list if only one url or scene given
    if type(urls) is str: urls = [urls]
    if type(scenes) is str: scenes = [scenes]

    ## run through urls and download the file
    if verbosity > 0: print('Downloading {} scenes'.format(len(urls)))
    lfiles, zfiles = [], []
    for ui, url in enumerate(urls):
        ## if scene list is provided use those scene names
        if len(scenes) == len(urls):
            scene = scenes[ui]
            session = None
        else:
            ## get access token - do it per url
            if verbosity > 1: print('Getting CDSE access token')
            data = {"client_id": "cdse-public", "grant_type":
                    "password","username": auth[0], "password": auth[1]}
            response = requests.post(auth_url, data=data, verify=True, allow_redirects=False)
            if 'access_token' not in response.json():
                print('Could not get access token for {}'.format(url))
                continue

            access_token = response.json()['access_token']
            if verbosity > 1: print(access_token)
            time.sleep(3)

            ## set up download session
            session = requests.Session()
            session.headers["Authorization"] = f"Bearer {access_token}"

            ## find out scene name
            response = requests.get(url.strip('/$value'))
            scene_atts = response.json()
            if 'Name' in scene_atts:
                scene = scene_atts['Name']
            if 'value' in scene_atts:
                scene = scene_atts['value'][0]['Name']

        ## local files
        lfile = '{}/{}'.format(output, scene)
        zfile = '{}/{}.zip'.format(output, scene[0:scene.find('.')])

        ## download if we don't have the scene
        if (override | (not os.path.exists(lfile)) & (not os.path.exists(zfile))):
            ## if scene list provided we need to set up the session here
            if session is None:
                if verbosity > 1: print('Getting CDSE access token')
                data = {"client_id": "cdse-public", "grant_type": "password","username": auth[0], "password": auth[1]}
                response = requests.post(auth_url, data=data, verify=True, allow_redirects=False)
                access_token = response.json()['access_token']
                if verbosity > 1: print(access_token)
                time.sleep(3)
                ## set up download session
                session = requests.Session()
                session.headers["Authorization"] = f"Bearer {access_token}"

            ## try url
            print('Downloading {}'.format(scene))
            print('Download URL: {}'.format(url))
            response = session.get(url, allow_redirects=False)
            ## follow redirects
            while response.status_code in (301, 302, 303, 307):
                url = response.headers['Location']
                response = session.get(url, allow_redirects=False)

            ## download file
            if os.path.exists(zfile): os.remove(zfile)
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
                if not os.path.exists(lfile):
                    print('Error extracting {}'.format(zfile))
                else:
                    print('Wrote {}'.format(lfile))

            ## remove downloaded zip file after extraction
            if remove_zip:
                if os.path.exists(zfile) & os.path.exists(lfile):
                    print('Deleting {}'.format(zfile))
                    os.remove(zfile)
                    print('Deleted {}'.format(zfile))

        ## list local paths
        if extract_zip:
            if os.path.exists(lfile): lfiles.append(lfile)
        else:
            if os.path.exists(zfile): lfiles.append(zfile)
    return(lfiles)
