## def azimuth_two_points
## computes azimuth between two points lon1,lat1 and lon2,lat2
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-04-14
## modifications: 2021-05-26 (QV) changed math to numpy

def azimuth_two_points(lon1, lat1, lon2, lat2):
    import numpy as np

    lat1r = np.radians(lat1)
    lat2r = np.radians(lat2)

    londiffr = np.radians(lon2-lon1)

    x = np.sin(londiffr) * np.cos(lat2r)
    y = np.cos(lat1r) * np.sin(lat2r) - (np.sin(lat1r)
                * np.cos(lat2r) * np.cos(londiffr))

    fw_azi = np.arctan2(x, y)
    fw_azi_deg = np.degrees(fw_azi)
    azimuth = (fw_azi_deg + 360.) % 360.

    return(azimuth)
