## def auth
## gets authentication token from USGS EarthExplorer API
## written by Quinten Vanhellemont, RBINS
## 2023-09-19
## modifications: 2023-09-25 (QV) changed to separate earthexplorer credentials

def auth(api_url = None, return_auth = False, netrc_machine = 'earthexplorer'):
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
        if ('EARTHEXPLORER_u' in os.environ) & \
           ('EARTHEXPLORER_p' in os.environ):
            auth = (os.environ['EARTHEXPLORER_u'], os.environ['EARTHEXPLORER_p'])

    ## get auth from config
    if auth is None:
        if (ac.config['EARTHEXPLORER_u'] != '') & \
           (ac.config['EARTHEXPLORER_p'] != ''):
            auth = (ac.config['EARTHEXPLORER_u'], ac.config['EARTHEXPLORER_p'])

    if auth is None:
        print('Could not determine EarthExplorer credentials. Please add them to your .netrc file or ACOLITE config file.')
        return

    ## get api URL from config
    if api_url is None: api_url = ac.config['EARTHEXPLORER_api']

    ## get access token
    data = {"username": auth[0], "password": auth[1]}
    response = requests.post(api_url+'/login', data = json.dumps(data))
    access_token = response.json()["data"]

    if return_auth: return(access_token, auth)

    return(access_token)
