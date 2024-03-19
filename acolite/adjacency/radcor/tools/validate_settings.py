## def validate_settings
## check if RAdCor settings make sense
##
## written by Alexandre Castagna, UGent (R version)
##            Quinten Vanhellemont, RBINS (Python version)
##
## for the RAdCor project 2023-2024
##
## modifications: 2024-03-18 (QV) added as function

def validate_settings(settings):
    ## initially assume settings are valid
    valid = True

    ## check forced aerosol model
    if settings['radcor_force_model'] is not None:
        modlist = ['M', 'C']
        if settings['radcor_force_model'] not in modlist:
            print('Error: forcemodel must be one of: {}'.format(', '.join(modlist)))
            valid = False

    ## check forced aerosol optical thickness
    if settings['radcor_force_aot'] is not None:
        if settings['radcor_force_aot'] < 0:
            print('Error: radcor_force_aot must not be negative.')
            valid = False
        if ~np.isfinite(settings['radcor_force_aot']):
            print('Error: radcor_force_aot must be finite.')
            valid = False

    ## check if edge extension method exists
    if settings['radcor_edge_extend']:
        expandlist = ['average', 'mirror', 'nearest']
        if settings['radcor_edge_extend_method'] not in expandlist:
            print('Error: radcor_edge_extend_method must be one of: {}'.format(', '.join(expandlist)))
            valid = False

    ## check how to complete the psf
    if not settings['radcor_psf_rescale']:
        completelist = ['average', 'neighborhood']
        if settings['radcor_psf_complete_method'] not in completelist:
            print('Error: radcor_psf_complete_method must be one of: {}'.format(', '.join(completelist)))
            valid = False

    return(valid)
