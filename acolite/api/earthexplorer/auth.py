## def auth
## gets authentication token from USGS EarthExplorer API
## written by Quinten Vanhellemont, RBINS
## 2023-09-19
## modifications: 2023-09-25 (QV) changed to separate earthexplorer credentials
##                2024-04-27 (QV) moved to acolite.api

def auth(api_url = None, return_auth = False, netrc_machine = 'earthexplorer'):
    import os, requests, json, netrc
    import acolite as ac

    ## get credentials
    auth = ac.shared.auth(netrc_machine)

    if auth is None:
        print('Could not determine {} credentials {}_u and {}_p.'.format(netrc_machine.upper(), netrc_machine.upper(), netrc_machine.upper()))
        print('Please add them to your .netrc file, environment variables, or ACOLITE config file.')
        return()

    ## get api URL from config
    if api_url is None: api_url = ac.config['EARTHEXPLORER_api']

    ## get access token
    data = {"username": auth[0], "password": auth[1]}
    response = requests.post(api_url+'/login', data = json.dumps(data))
    access_token = response.json()["data"]

    if return_auth: return(access_token, auth)

    return(access_token)
