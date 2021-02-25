## def geolocation_sub
## finds crop in lon and lat datasets
## written by Quinten Vanhellemont, RBINS 2021-02-25

def geolocation_sub(lat, lon, limit):
    import acolite as ac
    import numpy as np

    ## new version
    tmp = (lat >= limit[0]) & (lat <= limit[2]) & \
          (lon >= limit[1]) & (lon <= limit[3])
    roi = np.where(tmp)
    tmp = None

    if len(roi[0]) == 0:
        sub = None
    else:
        x0, x1 = min(roi[0]),max(roi[0])
        y0, y1 = min(roi[1]),max(roi[1])
        sub = [int(s) for s in [y0, x0, y1-y0, x1-x0]]

    return(sub)
