## def coef_band
## reads Nechad hyperspectral data and provides band specific coefficient
## either for band centre wavelength (set wave)
## written by Quinten Vanhellemont, RBINS
## 2023-07-30
## modifications:

def coef_band(nechad_par, cal_year = None, wave = None, sensor = None, band = None, rsr = None):
    import acolite as ac

    if nechad_par[0:3].upper() in ['T', 'TUR', 'TURBIDITY']:
        if cal_year is None: cal_year = 2009
    elif nechad_par[0:3].upper() in ['S', 'SPM', 'TSM']:
        if cal_year is None: cal_year = 2010
    else:
        print('Nechad parameter "{}" not configured.'.format(nechad_par))
        return

    ## load hyperspectral coefficient
    nechad_dict = ac.parameters.nechad.coef_hyper(nechad_par, year=cal_year)

    ## determine if specific wavelength is needed,
    ## or resampling to a specific band
    A_Nechad, C_Nechad = None, None
    if (sensor is None) & (band is None):
        if (wave is None):
            print('Please specify wavelength or sensor+band.')
        else:
            ## find closest wavelength in
            didx,algwave = ac.shared.closest_idx(nechad_dict['wave'], float(wave))
            A_Nechad = nechad_dict['A'][didx]
            C_Nechad = nechad_dict['C'][didx]
    else:
        if rsr is None:
            rsr = ac.shared.rsr_dict(sensor = sensor)[sensor]

        ## resample parameters to band
        #cdct = ac.shared.rsr_convolute_dict(nechad_dict['wave']/1000, nechad_dict['C'], rsrd[gem.gatts['sensor']]['rsr'])
        #adct = ac.shared.rsr_convolute_dict(nechad_dict['wave']/1000, 1/nechad_dict['A'], rsrd[gem.gatts['sensor']]['rsr'])
        #adct = {k:1/adct[k] for k in adct}
        #A_Nechad = adct[nechad_band]
        #C_Nechad = cdct[nechad_band]

        ## band specific rsr
        if 'rsr' in rsr.keys(): rsr = {band: rsr['rsr'][band]}

        ## resample parameters to band
        cdct = ac.shared.rsr_convolute_dict(nechad_dict['wave']/1000, nechad_dict['C'], rsr)
        adct = ac.shared.rsr_convolute_dict(nechad_dict['wave']/1000, 1/nechad_dict['A'], rsr)
        adct = {k:1/adct[k] for k in adct}
        A_Nechad = adct[band]
        C_Nechad = cdct[band]

    return(A_Nechad, C_Nechad)
