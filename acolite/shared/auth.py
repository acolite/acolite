## def auth
## gets authentication credentials
## written by Quinten Vanhellemont, RBINS
## 2023-09-25
## modifications: 2025-04-03 (QV) added run settings

def auth(machine):
    import os, requests, json, netrc
    import acolite as ac

    auth = None
    ## get auth from netrc file
    try:
        nr = netrc.netrc()
        ret = nr.authenticators(machine)
        if ret is not None:
            login, account, password = ret
            login = login.strip('"')
            password = password.strip('"')
            auth = (login, password)
    except:
        pass

    ## get auth from environment
    if auth is None:
        if ('{}_u'.format(machine.upper()) in os.environ) & \
           ('{}_p'.format(machine.upper()) in os.environ):
            auth = (os.environ['{}_u'.format(machine.upper())], \
                    os.environ['{}_p'.format(machine.upper())])

    ## get auth from config
    if auth is None:
        if ('{}_u'.format(machine.upper()) in ac.config) & \
           ('{}_p'.format(machine.upper()) in ac.config):
            if (ac.config['{}_u'.format(machine.upper())] != '') & \
               (ac.config['{}_p'.format(machine.upper())] != ''):
                auth = (ac.config['{}_u'.format(machine.upper())], \
                        ac.config['{}_p'.format(machine.upper())])

    ## get auth from run settings
    if auth is None:
        if ('{}_u'.format(machine.upper()) in ac.settings['run']) & \
           ('{}_p'.format(machine.upper()) in ac.settings['run']):
            if (ac.settings['run']['{}_u'.format(machine.upper())] != '') & \
               (ac.settings['run']['{}_p'.format(machine.upper())] != ''):
                auth = (ac.settings['run']['{}_u'.format(machine.upper())], \
                        ac.settings['run']['{}_p'.format(machine.upper())])

    if auth is None:
        print('Could not determine {} credentials. Please add them to your .netrc file or ACOLITE credentials file.'.format(machine))
        return

    return(auth)
