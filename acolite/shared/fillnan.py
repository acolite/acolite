## def fillnan
## fill array nans with closest values
## written by Quinten Vanhellemont
## 2020-10-28
## modifications: 2023-07-25 (QV) added max_distance keyword

def fillnan(data, max_distance=None):
    from scipy.ndimage import distance_transform_edt
    import numpy as np

    ## fill nans with closest value
    dis, ind = distance_transform_edt(np.isnan(data), return_distances=True, return_indices=True)
    data_filled = data[tuple(ind)]

    ## fill again with nan if greater than max_distance
    if (max_distance is not None):
        if max_distance > 0:
            sub = np.where(dis > max_distance)
            data_filled[sub] = np.nan

    return(data_filled)
