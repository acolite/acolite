## def validate_settings
## check if RAdCor settings make sense
##
## written by Alexandre Castagna, UGent (R version)
##            Quinten Vanhellemont, RBINS (Python version)
##
## for the RAdCor project 2023-2024
##
## modifications: 2024-03-18 (QV) added as function
##                2024-03-25 (QV) added numpy import
##                2024-03-26 (QV) added check for tsdsf_psf_complete_method

def validate_settings(settings):
    import numpy as np

    ## initially assume settings are valid
    valid = True

    ## check aot estimate
    estimatelist = ['tsdsf', 'dsf', 'optimise']
    if settings['radcor_aot_estimate'] not in estimatelist:
        print('Error: radcor_aot_estimate must be one of: {}'.format(', '.join(estimatelist)))
        valid = False
    if settings['radcor_aot_estimate'] == 'optimise':
        for k in ['radcor_optimise_target_lon', 'radcor_optimise_target_lat', 'radcor_optimise_target_rhos']:
            if (settings[k] is None):
                print('Provide {} for radcor_aot_estimate=optimise'.format(k))
                valid = False

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
        completelist = ['average', 'neighbourhood']
        if settings['radcor_psf_complete_method'] not in completelist:
            print('Error: radcor_psf_complete_method must be one of: {}'.format(', '.join(completelist)))
            valid = False

    ## check how to complete the psf
    if not settings['tsdsf_psf_rescale']:
        completelist = ['average', 'neighbourhood']
        if settings['tsdsf_psf_complete_method'] not in completelist:
            print('Error: tsdsf_psf_complete_method must be one of: {}'.format(', '.join(completelist)))
            valid = False

    return(valid)
