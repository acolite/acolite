## def auth
## gets authentication credentials
## written by Quinten Vanhellemont, RBINS
## 2023-09-25
## modifications: 2025-04-03 (QV) added run settings
##                2025-07-09 (QV) added token option


def auth(machine):
    import netrc
    import os

    import acolite as ac

    auth = None
    ## get auth from netrc file or NETRC environment variable
    try:
        nr = netrc.netrc(os.environ['NETRC'])
    except KeyError:
        nr = netrc.netrc()
    try:
        ret = nr.authenticators(machine)
        if ret is not None:
            login, account, password = ret
            login = login.strip('"')
            password = password.strip('"')
            auth = (login, password)
    except:
        pass

    ## remove token from machine if getting from environment/config/settings
    token = False
    if machine.lower().endswith('_token'):
        machine = machine[0:-6]
        token = True

    ## get auth from environment
    if auth is None:
        if ('{}_u'.format(machine.upper()) in os.environ):
            if ('{}_p'.format(machine.upper()) in os.environ):
                auth = (os.environ['{}_u'.format(machine.upper())], \
                        os.environ['{}_p'.format(machine.upper())])
            if (token) & ('{}_token'.format(machine.upper()) in os.environ):
                auth = (os.environ['{}_u'.format(machine.upper())], \
                        os.environ['{}_token'.format(machine.upper())])

    ## get auth from config
    if auth is None:
        if ('{}_u'.format(machine.upper()) in ac.config):
            if ('{}_p'.format(machine.upper()) in ac.config):
                if (ac.config['{}_u'.format(machine.upper())] != '') & \
                   (ac.config['{}_p'.format(machine.upper())] != ''):
                    auth = (ac.config['{}_u'.format(machine.upper())], \
                            ac.config['{}_p'.format(machine.upper())])
            if (token) & ('{}_token'.format(machine.upper()) in ac.config):
                if (ac.config['{}_u'.format(machine.upper())] != '') & \
                   (ac.config['{}_token'.format(machine.upper())] != ''):
                    auth = (ac.config['{}_u'.format(machine.upper())], \
                            ac.config['{}_token'.format(machine.upper())])

    ## get auth from run settings
    if auth is None:
        if ('{}_u'.format(machine.upper()) in ac.settings['run']):
            if ('{}_p'.format(machine.upper()) in ac.settings['run']):
                if (ac.settings['run']['{}_u'.format(machine.upper())] != '') & \
                   (ac.settings['run']['{}_p'.format(machine.upper())] != ''):
                    auth = (ac.settings['run']['{}_u'.format(machine.upper())], \
                            ac.settings['run']['{}_p'.format(machine.upper())])
            if (token) & ('{}_token'.format(machine.upper()) in ac.settings['run']):
                if (ac.settings['run']['{}_u'.format(machine.upper())] != '') & \
                   (ac.settings['run']['{}_token'.format(machine.upper())] != ''):
                    auth = (ac.settings['run']['{}_u'.format(machine.upper())], \
                            ac.settings['run']['{}_token'.format(machine.upper())])

    if auth is None:
        print('Could not determine {} credentials. Please add them to your .netrc file or ACOLITE credentials file.'.format(machine))
        return

    return(auth)
