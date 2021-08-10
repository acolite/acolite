## Read DESIS xml metadata file
## written by Quinten Vanhellemont, RBINS
## 2021-08-10

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
    tags = ['startTime','endTime', 'mission', 'satelliteID', 'sensor', 'productType']
    tags += ['sunAzimuthAngle', 'sunZenithAngle', 'sceneAzimuthAngle', 'sceneIncidenceAngle']
    tags += ['pixelSize', 'widthOfScene', 'heightOfScene']

    for tag in tags:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) > 0: meta[tag] = tdom[0].firstChild.nodeValue

    return(meta)
