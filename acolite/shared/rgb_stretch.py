## def rgb_stretch
## does stretching for RGB visualisation
## written by Quinten Vanhellemont, RBINS
## 2022-12-30
## modifications: 2022-12-30 (QV) moved gamma stretching

def rgb_stretch(data, gamma=1.0, stretch='linear', bsc=None, percentiles=(5,95)):
    import numpy as np

    stretches = ['linear', 'log', 'sinh', 'sqrt']
    if stretch not in stretches:
        print('rgb_stretch={} not configured'.format(setu['rgb_stretch']))
        print('Using default linear stretch')
        stretch = 'linear'

    ## set scale range
    if bsc is None: bsc = np.nanpercentile(data, percentiles)

    ## do stretch
    data_ = data**gamma
    if stretch == 'linear':
        data_ = np.interp(data_, bsc, [0, 1])
    elif stretch == 'log':
        if bsc[0] <= 0: bsc[0] = 0.01
        data_ = np.interp(np.log(data_), np.log(bsc) , [0, 1])
    elif stretch == 'sinh':
        data_ = np.interp(np.sinh(data_), np.sinh(bsc) , [0, 1])
    elif stretch == 'sqrt':
        data_ = np.interp(np.sqrt(data_), np.sqrt(bsc) , [0, 1])
    return(data_)
