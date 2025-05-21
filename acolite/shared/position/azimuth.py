## def azimuth
## compute azimuth between two lon, lat positions
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-20
## modifications:

def azimuth(lon1, lat1, lon2, lat2):
    import numpy as np

    lon1 = np.radians(lon1)
    lon2 = np.radians(lon2)
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)

    x = np.sin(lon2-lon1) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(lon2-lon1)
    azi = np.degrees(np.arctan2(x, y))
    azi = np.atleast_1d(azi)

    ## shift negatives
    azi[azi<0] += 360

    ## return single value if one position was asked
    if azi.shape == (1,): azi = azi[0]
    return(azi)
