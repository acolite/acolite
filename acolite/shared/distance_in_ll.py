## approximates the distance in one degree at the given latitude
##
## written by Quinten Vanhellemont
## 2018-02-07
## modifications:

def distance_in_ll(lat=0, earth_radius=6378.1370):
    from numpy import pi, cos
    circle_lon=2.*pi*earth_radius # length of one meridian
    onedeglat=(circle_lon/360.) # km in one degree of latitude

    proj_earth_radius=earth_radius*cos(lat*pi/180.)
    circle_lat=2.*pi*proj_earth_radius # length of circle of mean longitude
    onedeglon=(circle_lat/360.) # km in one degree of longitude

    return(onedeglon, onedeglat)
