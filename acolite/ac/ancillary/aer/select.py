## select
## select model and aot based on ancillary aerosol data
##
## written by Quinten Vanhellemont, RBINS
## 2025-12-08
## modifications:

def select(date, lon, lat):

    import acolite as ac
    import numpy as np

    ## retrieve ancillary aerosols
    aer_anc = ac.ac.ancillary.aer.get(date, lon, lat)

    if 'data' not in aer_anc:
        print('Error in ancillary AOT retrieval.')
        return
    elif aer_anc['data'] == {}:
        print('Error in ancillary AOT retrieval.')
        return
    else:
        aer_ang = aer_anc['data']['TOTANGSTR']['interp']
        aer_aot = aer_anc['data']['TOTEXTTAU']['interp']
        aer_ang_mean = np.nanmean(aer_ang)

        ## 6SV models C(ontinental) 1.08, M(aritime) 0.28,
        anc_ang_threshold = 0.68 ## midway between C and M
        anc_lut = ac.settings['run']['luts'][0]
        if len(ac.settings['run']['luts']) > 1:
            print('Selecting LUT based on ancillary mean angstrom {:.2f} from LUTs {}'.format(aer_ang_mean, ac.settings['run']['luts']))
            print('Angstrom threshold M < {} < C'.format(anc_ang_threshold))
            for lut in ac.settings['run']['luts']:
                if (aer_ang_mean >= anc_ang_threshold) & (lut[-1] == '1'): anc_lut = '{}'.format(lut)
                elif (aer_ang_mean < anc_ang_threshold) & (lut[-1] == '2'): anc_lut = '{}'.format(lut)
        anc_aot = np.nanmean(aer_aot) * 1
    print('Setting dsf_fixed_aot={:.3f} (mean) and dsf_fixed_lut={} (mean angstrom={:.2f}) based on ancillary data'.format(anc_aot, anc_lut, aer_ang_mean))

    return(anc_lut, anc_aot)
