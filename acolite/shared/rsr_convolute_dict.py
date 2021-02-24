## def rsr_convolute_dict
## resample given data to sensor rsr
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-11
## modifications: 2017-10-17 (QV) moved wave range and step to keywords, renamed from interprsr
##                2020-03-02 (QV) added ceil to number of elements in linspace

def rsr_convolute_dict(wave_data, data, rsr, wave_range=[0.2,2.4], wave_step=0.001):
    import numpy as np

    ## set up wavelength space
    wave_hyper = np.linspace(wave_range[0],wave_range[1],int(((wave_range[1]-wave_range[0])/wave_step)+2))

    ## interpolate RSR to same dimensions
    rsr_hyper = dict()
    for band in rsr:
        band_wave_hyper = wave_hyper
        band_response_hyper = np.interp(wave_hyper, rsr[band]['wave'], rsr[band]['response'], left=0, right=0)
        band_response_sum = sum(band_response_hyper)
        rsr_hyper[band]={'wave':band_wave_hyper, 'response': band_response_hyper, 'sum':band_response_sum}

    resdata={}
    data_hyper = np.interp(wave_hyper, wave_data, data, left=0, right=0)
    for band in rsr:
        resdata[band] = (sum(data_hyper*rsr_hyper[band]['response'])/rsr_hyper[band]['sum'])
    return resdata
