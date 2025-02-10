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
##                2024-03-26 (QV) added check for tsdsf_kernel_complete_method
##                2024-05-21 (QV) added radcor_optimise_aot_cost and radcor_optimise_target_type
##                2024-12-16 (QV) removed radcor/tsdsf_kernel_rescale and added renormalise to radcor/tsdsf_kernel_complete_method
##                2025-02-10 (QV) renamed radcor_optimise_ to optimise_

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
        for k in ['optimise_target_lon', 'optimise_target_lat', 'optimise_aot_cost']:
            if (settings[k] is None):
                print('Provide {} for radcor_aot_estimate=optimise'.format(k))
                valid = False
            if (settings['optimise_target_rhos'] is None) & (settings['optimise_target_rhos_file'] is None):
                print('Provide either optimise_target_rhos or optimise_target_rhos_file for radcor_aot_estimate=optimise')
                valid = False
        costlist = ['RMSD', 'MAPD']
        if (valid) & (settings['optimise_aot_cost'].upper() not in costlist):
            print('Error: optimise_aot_cost must be one of: {}'.format(', '.join(costlist)))
            valid = False
        targetlist = ['pixel', 'box', 'circle']
        if (valid) & (settings['optimise_target_type'] not in targetlist):
            print('Error: optimise_target_type must be one of: {}'.format(', '.join(targetlist)))
            valid = False
        if (valid) & (settings['optimise_target_type'] != 'pixel') & (settings['optimise_target_units'][0].lower() not in ['p', 'm']):
             print('Error: optimise_target_units must be m(etre) or p(ixel)')
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

    ## check how to complete the psf for radcor
    completelist = ['average', 'neighbourhood', 'renormalise']
    if settings['radcor_kernel_complete_method'] not in completelist:
        print('Error: radcor_kernel_complete_method must be one of: {}'.format(', '.join(completelist)))
        valid = False

    ## check how to complete the psf for tsdsf
    completelist = ['average', 'neighbourhood', 'renormalise']
    if settings['tsdsf_kernel_complete_method'] not in completelist:
        print('Error: tsdsf_kernel_complete_method must be one of: {}'.format(', '.join(completelist)))
        valid = False

    return(valid)
