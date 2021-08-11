## def metadata
## read GF6 metadata
##
## written by Quinten Vanhellemont, RBINS
## 2021-07-30
## modifications: 2021-08-09 (QV) function added to generic acolite

def metadata(metafile):
    import dateutil.parser
    from xml.dom import minidom

    import acolite as ac
    import numpy as np

    try:
        xmldoc = minidom.parse(metafile)
    except:
        print('Error opening metadata file.')
        sys.exit()

    meta = {}
    tags = ['SatelliteID','SensorID', 'ProductLevel', 'StartTime', 'EndTime', 'CenterTime', 'Bands']
    tags += ['SolarAzimuth', 'SolarZenith', 'SatelliteAzimuth', 'SatelliteZenith']
    tags += ['WidthInPixels', 'HeightInPixels', 'WidthInMeters', 'HeightInMeters']

    tags += ['CenterLatitude', 'CenterLongitude',
            'TopLeftLatitude', 'TopLeftLongitude',
            'TopRightLatitude', 'TopRightLongitude',
            'BottomRightLatitude', 'BottomRightLongitude',
            'BottomLeftLatitude', 'BottomLeftLongitude']

    band_names = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8' ]
    band_names += ['MS1', 'MS2', 'MS3', 'MS4', 'PAN']

    tags += band_names

    for tag in tags:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) > 0: meta[tag] = tdom[0].firstChild.nodeValue

    return(meta)
