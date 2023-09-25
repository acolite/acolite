## def auth
## gets authentication token from USGS EarthExplorer API
## written by Quinten Vanhellemont, RBINS
## 2023-09-19
## modifications:

def auth(api_url = None, return_auth = False, netrc_machine = 'earthdata'):
    import os, requests, json, netrc
    import acolite as ac

    auth = None

    ## get auth from netrc file
    try:
        nr = netrc.netrc()
        ret = nr.authenticators(netrc_machine)
        if ret is not None:
            login, account, password = ret
            login = login.strip('"')
            password = password.strip('"')
            auth = (login, password)
    except:
        pass

    ## get auth from environment
    if auth is None:
        if ('EARTHDATA_u' in os.environ) & \
           ('EARTHDATA_p' in os.environ):
            auth = (os.environ['EARTHDATA_u'], os.environ['EARTHDATA_p'])

    ## get auth from config
    if auth is None:
        if (ac.config['EARTHDATA_u'] != '') & \
           (ac.config['EARTHDATA_p'] != ''):
            auth = (ac.config['EARTHDATA_u'], ac.config['EARTHDATA_p'])

    if auth is None:
        print('Could not determine EarthData credentials. Please add them to your .netrc file or ACOLITE config file.')
        return

    ## get api URL from config
    if api_url is None: api_url = ac.config['EARTHEXPLORER_api']

    ## get access token
    data = {"username": auth[0], "password": auth[1]}
    response = requests.post(api_url+'/login', data = json.dumps(data))
    access_token = response.json()["data"]

    if return_auth: return(access_token, auth)

    return(access_token)
