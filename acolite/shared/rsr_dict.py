## def rsr_dict
## imports rsr file(s)
## written by Quinten Vanhellemont, RBINS
## 2021-02-27
## modifications: 2022-02-25 (QV) added warning printout

def rsr_dict(sensor = None, rsrd = None):
    import glob, os
    import numpy as np
    import acolite as ac

    ## setup wavelength array
    waves = np.arange(250, 2502.5)/1000

    ## find rsr files
    if rsrd is None:
        if sensor is None:
            sens = glob.glob(ac.config['data_dir']+'/RSR/*.txt')
        else:
            sens = glob.glob(ac.config['data_dir']+'/RSR/{}.txt'.format(sensor))
            if len(sens) == 0: print('Could not find {} RSR file at {}'.format(sensor, ac.config['data_dir']+'/RSR/'))

        rsrd = {}
        for rsrf in sens:
            ## read rsr file
            fsensor = os.path.basename(rsrf).split('.txt')[0]
            rsr, rsr_bands = ac.shared.rsr_read(rsrf)
            rsrd[fsensor] = {'rsr':rsr, 'rsr_bands':rsr_bands}

    for fsensor in rsrd:
        ## compute band weighted wavelengths and band names
        rsrd[fsensor]['wave_mu'] = ac.shared.rsr_convolute_dict(waves, waves, rsrd[fsensor]['rsr'])
        rsrd[fsensor]['wave_nm'] = {b:rsrd[fsensor]['wave_mu'][b]*1000 for b in rsrd[fsensor]['wave_mu']}
        rsrd[fsensor]['wave_name'] = {b:'{:.0f}'.format(rsrd[fsensor]['wave_nm'][b]) for b in rsrd[fsensor]['wave_nm']}
        if 'rsr_bands' not in rsrd[fsensor]:
            rsrd[fsensor]['rsr_bands'] = [b for b in rsrd[fsensor]['rsr']]
    return(rsrd)
