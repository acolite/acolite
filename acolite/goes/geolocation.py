## def geolocation
## compute ABI geolocation based on imager projection info in NetCDF file
##
## written by Quinten Vanhellemont, RBINS
## 2025-06-05
## modifications:

def geolocation(file):
    import acolite as ac
    import numpy as np

    ## open file
    gemi = ac.gem.gem(file)
    ## get projection attributes
    goes_imager_projection = gemi.atts('goes_imager_projection')
    ## get scanning and elevation angles
    x = gemi.data('x')
    y = gemi.data('y')
    gemi.close()
    gemi = None

    ## create 2d x and y coordinates
    x2d, y2d = np.meshgrid(x, y)

    ## satellite heigth
    height = goes_imager_projection['perspective_point_height'] + goes_imager_projection['semi_major_axis']

    ## satellite longitude
    lambda_0 = np.radians(goes_imager_projection['longitude_of_projection_origin'])

    #a = np.power(np.sin(x2d),2.0) + (np.power(np.cos(x2d),2.0) * (np.power(np.cos(y2d),2.0) + (((goes_imager_projection['semi_major_axis'] ** 2)/(goes_imager_projection['semi_minor_axis'] ** 2))*np.power(np.sin(y2d),2.0))))
    a = np.sin(x2d) ** 2.0 + (np.cos(x2d) ** 2.0 * np.cos(y2d) ** 2.0 + (((goes_imager_projection['semi_major_axis'] ** 2) / (goes_imager_projection['semi_minor_axis'] ** 2)) * np.sin(y2d) ** 2.0))
    b = -2.0 * height * np.cos(x2d) * np.cos(y2d)
    c = (height ** 2.0) - (goes_imager_projection['semi_major_axis'] ** 2.0)

    ## compute satellite position
    rs = (-1.0 * b - np.sqrt((b ** 2) - (4.0 * a * c)))/(2.0 * a)
    sx = rs * np.cos(x2d) * np.cos(y2d)
    sy = -rs * np.sin(x2d)
    sz = rs * np.cos(x2d) * np.sin(y2d)
    del rs, x2d, y2d

    ## compute lon and lat
    lon = np.degrees(lambda_0 - np.arctan(sy / (height - sx)))
    lat = np.degrees(np.arctan(((goes_imager_projection['semi_major_axis'] ** 2)/(goes_imager_projection['semi_minor_axis'] ** 2)) * ((sz / np.sqrt(((height - sx) ** 2) + (sy ** 2))))))
    del sx, sy, sz

    return(lon, lat)
