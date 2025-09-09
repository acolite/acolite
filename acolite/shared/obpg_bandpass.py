## def obpg_bandpass
## get rsr bandpass file from OBPG
## written by Quinten Vanhellemont, RBINS
## 2025-09-09
## modifications:

def obpg_bandpass(sensor, local_dir = None, url = 'https://oceancolor.gsfc.nasa.gov/images/rsr/'):
    import os
    import numpy as np
    import acolite as ac
    if sensor not in ['oci', 'oci_l1b']:
         print('Sensor {} to be configured.'.format(sensor))
         return

    if local_dir is None: local_dir = '{}/OBPG'.format(ac.config['external_dir'])

    file = '{}_bandpass.csv'.format(sensor)
    local_file = '{}/{}'.format(local_dir, file)

    if not os.path.exists(local_file):
        remote_file = '{}/{}'.format(url, file)
        print('Downloading {} to {}'.format(remote_file, local_file))
        ac.shared.download_file(remote_file, local_file)

    ## simple import for now
    if os.path.exists(local_file):
        rsr_bandpass = np.loadtxt(local_file, skiprows = 2, delimiter = ',')
        return(rsr_bandpass)
