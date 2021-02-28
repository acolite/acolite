## def rsr_dict
## imports rsr file(s)
## written by Quinten Vanhellemont, RBINS
## 2021-02-27
## modifications:

def rsr_dict(sensor=None):
    import glob, os
    import numpy as np
    import acolite as ac

    ## setup wavelength array
    waves = np.arange(250, 2502.5)/1000

    ## find rsr files
    if sensor is None:
        sens = glob.glob(ac.config['data_dir']+'/RSR/*.txt')
    else:
        sens = glob.glob(ac.config['data_dir']+'/RSR/{}.txt'.format(sensor))

    rsrd = {}
    for rsrf in sens:
        ## read rsr file
        fsensor = os.path.basename(rsrf).split('.txt')[0]
        rsr, rsr_bands = ac.shared.rsr_read(rsrf)
        rsrd[fsensor] = {'rsr':rsr, 'rsr_bands':rsr_bands}

        ## compute band weighted wavelengths and band names
        rsrd[fsensor]['waves_mu'] = ac.shared.rsr_convolute_dict(waves, waves, rsrd[fsensor]['rsr'])
        rsrd[fsensor]['waves_nm'] = {b:rsrd[fsensor]['waves_mu'][b]*1000 for b in rsrd[fsensor]['waves_mu']}
        rsrd[fsensor]['waves_name'] = {b:'{:.0f}'.format(rsrd[fsensor]['waves_nm'][b]) for b in rsrd[fsensor]['waves_nm']}

    return(rsrd)
