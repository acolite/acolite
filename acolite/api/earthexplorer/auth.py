## def auth
## gets authentication token from USGS EarthExplorer API
## written by Quinten Vanhellemont, RBINS
## 2023-09-19
## modifications: 2023-09-25 (QV) changed to separate earthexplorer credentials
##                2024-04-27 (QV) moved to acolite.api
##                2025-04-03 (QV) use login-token

def auth(api_url = None, return_auth = False, netrc_machine = 'earthexplorer',
         user = None, password = None, token = None):
    import os, requests, json, netrc
    import acolite as ac

    ## get credentials
    auth = ac.shared.auth(netrc_machine)

    if auth is None:
        print('Could not determine {} credentials {}_u and {}_p.'.format(netrc_machine.upper(), netrc_machine.upper(), netrc_machine.upper()))
        print('Please add them to your .netrc file, environment variables, or ACOLITE config file.')

    ## get username
    # user_key = 'EARTHEXPLORER_u'
    # if (user is None) & (user_key in ac.settings['run']):
    #     if (ac.settings['run'][user_key] not in [None, '']):
    #         user = '{}'.format(ac.settings['run'][user_key])
    # if (user is None) & (user_key in ac.config):
    #     if (ac.config[user_key] not in [None, '']):
    #         user = '{}'.format(ac.config[user_key])
    # if (user is None) & (user_key in os.environ):
    #     if (os.environ[user_key] not in [None, '']):
    #         user = '{}'.format(os.environ[user_key])
    # if (user is None):
    #     print('Could not determine {}.'.format(user_key))
    #     print('Please provide an {} as environment variable, or in the settings or credentials file.'.format(user_key))
    #     return()

    ## get password
    # pass_key = 'EARTHEXPLORER_p'
    # if (password is None) & (pass_key in ac.settings['run']):
    #     if (ac.settings['run'][pass_key] not in [None, '']):
    #         password = '{}'.format(ac.settings['run'][pass_key])
    # if (password is None) & (pass_key in ac.config):
    #     if (ac.config[pass_key] not in [None, '']):
    #         password = '{}'.format(ac.config[pass_key])
    # if (password is None) & (pass_key in os.environ):
    #     if (os.environ[pass_key] not in [None, '']):
    #         password = '{}'.format(os.environ[pass_key])
    # if (password is None):
    #     print('Could not determine {}.'.format(pass_key))
    #     print('Please provide an {} as environment variable, or in the settings or credentials file.'.format(pass_key))
    #     return()

    ## get token
    token_key = 'EARTHEXPLORER_token'
    if (token is None) & (token_key in ac.settings['run']):
        if (ac.settings['run'][token_key] not in [None, '']):
            token = '{}'.format(ac.settings['run'][token_key])
    if (token is None) & (token_key in ac.config):
        if (ac.config[token_key] not in [None, '']):
            token = '{}'.format(ac.config[token_key])
    if (token is None) & (token_key in os.environ):
        if (os.environ[token_key] not in [None, '']):
            token = '{}'.format(os.environ[token_key])
    if (token is None):
        print('Could not determine {}.'.format(token_key))
        print('Please provide an {} as environment variable, or in your settings or the credentials file.'.format(token_key))
        return()

    ## get api URL from config
    if api_url is None: api_url = ac.config['EARTHEXPLORER_api']

    ## get access token
    data = {"username": auth[0], "token": token}
    response = requests.post(api_url+'/login-token', data = json.dumps(data))
    access_token = response.json()["data"]

    if return_auth: return(access_token, (auth[0], auth[1], token))

    return(access_token)
