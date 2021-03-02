## def azimuth_two_points
## computes azimuth between two points lon1,lat1 and lon2,lat2
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-04-14
## modifications:

def azimuth_two_points(lon1, lat1, lon2, lat2):
    from math import radians, pow, sin, cos, atan2, degrees
    lat1r = radians(lat1)
    lat2r = radians(lat2)

    londiffr = radians(lon2-lon1)

    x = sin(londiffr) * cos(lat2r)
    y = cos(lat1r) * sin(lat2r) - (sin(lat1r)
                * cos(lat2r) * cos(londiffr))

    fw_azi = atan2(x, y)
    fw_azi_deg = degrees(fw_azi)
    azimuth = (fw_azi_deg + 360.) % 360.

    return azimuth
