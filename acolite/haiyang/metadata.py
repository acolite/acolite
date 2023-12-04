## def metadata
## read Haiyang CZI metadata
##
## written by Quinten Vanhellemont, RBINS
## 2023-02-19
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

    meta = {}
    tags = ['SatelliteID', 'SensorID', 'ProductLevel',
            'StartTime', 'EndTime', 'CentreTime', 'Bands']

    for tag in tags:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) > 0:
            if tdom[0].firstChild is None: continue
            meta[tag] = tdom[0].firstChild.nodeValue

    mtag = 'CentreLocation'
    cdom =  xmldoc.getElementsByTagName(mtag)
    for tag in ['Latitude', 'Longitude']:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) > 0:
            if tdom[0].firstChild is None: continue
            meta['{}_{}'.format(mtag, tag)] = float(tdom[0].firstChild.nodeValue)
    return(meta)
