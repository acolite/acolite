## def fillnan
## fill array nans with closest values
## written by Quinten Vanhellemont
## 2020-10-28
## modifications:

def fillnan(data):
    from scipy.ndimage import distance_transform_edt
    import numpy as np

    ## fill nans with closest value
    ind = distance_transform_edt(np.isnan(data), return_distances=False, return_indices=True)
    return(data[tuple(ind)])
