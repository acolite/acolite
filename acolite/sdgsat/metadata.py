## def metadata
## read SDBSAT1 KX10 metadata
##
## written by Quinten Vanhellemont, RBINS
## 2023-01-03
## modifications:

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

    metadata = {}
    tags = ['SatelliteID', 'SensorID',
            'RollSatelliteAngle','PitchSatelliteAngle', 'YawSatelliteAngle',
            'SolarAzimuth', 'SolarZenith']
    for tag in tags:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) > 0:
            if tdom[0].firstChild is None: continue
            metadata[tag] = tdom[0].firstChild.nodeValue

    timetags = ['StartTime', 'CenterTime', 'EndTime']
    for tt in timetags:
        cdom = xmldoc.getElementsByTagName(tt)
        if len(cdom) > 0:
            for cam in ['Acamera', 'Bcamera']:
                tag = '{}-{}'.format(tt, cam)
                tdom = cdom[0].getElementsByTagName(cam)
                if len(tdom) > 0:
                    if tdom[0].firstChild is None: continue
                    metadata[tag] = tdom[0].firstChild.nodeValue

    return(metadata)
