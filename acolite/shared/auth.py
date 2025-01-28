## def auth
## gets authentication credentials
## written by Quinten Vanhellemont, RBINS
## 2023-09-25
## modifications:

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
            auth = (os.environ['{}_u'.format(machine.upper())], os.environ['{}_p'.format(machine.upper())])

    ## get auth from config
    machine_to_table = {
        "earhtdata": "EarthData",
        "earthexplorer": "EarthExplorer",
        "cdse": "CDSE",
    }
    if auth is None:
        if ac.config[machine_to_table[machine]]["u"] and ac.config[machine_to_table[machine]]["p"]:
            auth = (ac.config[machine_to_table[machine]].values())

    if auth is None:
        print('Could not determine {} credentials. Please add them to your .netrc file or ACOLITE config file.'.format(machine))
        return

    return(auth)
