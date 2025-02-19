## def rsr_convolute_dict
## resample given data to sensor rsr
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-11
## modifications: 2017-10-17 (QV) moved wave range and step to keywords, renamed from interprsr
##                2020-03-02 (QV) added ceil to number of elements in linspace
##                2022-10-25 (QV) changed upper bound of default range to 2.5
##                2024-11-04 (QV) changed to np.nansum
##                2025-02-10 (QV) added fill_value keyword, check if nans

def rsr_convolute_dict(wave_data, data, rsr, wave_range=[0.2,2.55], wave_step=0.001, fill_value = 0):
    import numpy as np

    ## set up wavelength space
    wave_hyper = np.linspace(wave_range[0],wave_range[1],int(((wave_range[1]-wave_range[0])/wave_step)+1))

    ## interpolate RSR to same dimensions
    rsr_hyper = dict()
    for band in rsr:
        band_wave_hyper = wave_hyper
        band_response_hyper = np.interp(wave_hyper, rsr[band]['wave'], rsr[band]['response'], left=fill_value, right=fill_value)
        band_response_sum = np.nansum(band_response_hyper)
        rsr_hyper[band]={'wave':band_wave_hyper, 'response': band_response_hyper, 'sum':band_response_sum}

    resdata={}
    data_hyper = np.interp(wave_hyper, wave_data, data, left=fill_value, right=fill_value)
    for band in rsr:
        d = data_hyper*rsr_hyper[band]['response']
        if all(np.isnan(d)):
            resdata[band] = np.nan
        else:
            resdata[band] = (np.nansum(d)/rsr_hyper[band]['sum'])
    return resdata
