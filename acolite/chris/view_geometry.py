### chris_oaa_ova
## code from Héloïse Lavigne, RBINS
## adapted by QV on 2018-03-22
##
## modifications 2020-05-20 (QV) new version from HL
##               2021-06-08 (QV) renamed from chris_oaa_ova2 integrated in generic acolite

#### Python function to estimate Observation Zenith Angle and Observation Azimuth Angles when they are not available

# mza : minimum zenith Angle (degrees)
# flyby_za : Fly-by Zenith Angle (degrees)
# Palt : Platform Altitude (km)
# Talt : Target Altitude (km)
# Tlat : Target Latitude (degrees)

# return [Observation Zenith Angle, Observation Azimuth Angle] (degrees)
# code produced from the excel sheet CHRISProba-Obs-Angles-v2 created by H. Peter White on March 2016

def view_geometry(mza, flyby_za, Palt, Talt, Tlat, i=97.7162):
    import numpy as np

    flyby_za = -flyby_za
    i = i*np.pi/180 #rad
    Rearth = 6371.0088 # km
    mzar = mza*np.pi/180
    flyby_zar = flyby_za*np.pi/180
    Tlatr = Tlat*np.pi/180

    A = flyby_zar - np.arcsin((Rearth/(Rearth+Palt))*np.sin(flyby_zar))
    B = mzar - np.arcsin( (Rearth+Talt)/(Rearth+Palt)*np.sin(mzar))
    C = np.arccos(np.cos(A)*np.cos(B))
    D = np.arctan(np.sin(C)/(((Rearth+Palt)/Rearth) - np.cos(C)) )
    E = np.arctan(np.tan(i-np.pi/2)/np.sqrt(1 - (np.sin(Tlatr)/np.cos(i-np.pi/2))**2))

    RAAr = 0

    if C != 0 and  np.sin(A)/np.sin(C) <= -1 :
        RAAr = -np.pi/2

    if C != 0 and -1 < np.sin(A)/np.sin(C) and np.sin(A)/np.sin(C) <= 1 :
        RAAr = np.arcsin(np.sin(A)/np.sin(C))

    if C != 0 and 1 < np.sin(A)/np.sin(C) :
        RAAr = np.pi/2

    OA = E + np.pi/2

    if mzar < 0 :
        OA  = E - np.pi/2

    OZAr = C + D
    OAAr = OA - RAAr

    if mzar < 0:
        OAAr = 2*np.pi + OA + RAAr

    OZA = OZAr*180/np.pi
    OAA = OAAr*180/np.pi
    OAA = 360 - OAA

    return(OZA, OAA)
