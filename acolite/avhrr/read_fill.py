## def read_fill
## parses reads and fills gaps for AVHRR L1B/1C European Data Set bundle
## to do: add other interpolation methods
## written by Quinten Vanhellemont, RBINS
## 2024-11-14
## modifications:

def read_fill(image_file, ds, sub = None, method = 'nearest'):
    import numpy as np
    import acolite as ac

    data = ac.shared.nc_data(image_file, ds, sub = sub)
    data[data.mask] = np.nan
    if method == 'nearest':
        data = ac.shared.fillnan(data)

    return(data)
