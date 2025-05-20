## def metadata
## read DEIMOS2/GEOSAT2 metadata
##
## written by Quinten Vanhellemont, RBINS
## 2024-05-16
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
    tags = ['DATASET_NAME', 'GEOMETRIC_PROCESSING', 'SOURCE_ID',
            'NCOLS', 'NROWS', 'NBANDS',
            'MISSION', 'MISSION_INDEX', 'INSTRUMENT',
            'IMAGING_DATE', 'IMAGING_TIME', 'START_TIME', 'STOP_TIME', 'EARTH_SUN_DISTANCE',
            'VIEWING_ANGLE', 'INCIDENCE_ANGLE', 'SUN_ELEVATION', 'SUN_AZIMUTH',
            'PIXEL_RESOLUTION_X', 'PIXEL_RESOLUTION_Y']

    for tag in tags:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) > 0:
            if tdom[0].firstChild is None: continue
            meta[tag] = tdom[0].firstChild.nodeValue
            if tag in ['NCOLS', 'NROWS', 'NBANDS']:
                meta[tag] = int(meta[tag])
            if tag in ['EARTH_SUN_DISTANCE', 'VIEWING_ANGLE', 'INCIDENCE_ANGLE',
                       'SUN_ELEVATION', 'SUN_AZIMUTH', 'PIXEL_RESOLUTION_X', 'PIXEL_RESOLUTION_Y']:
                meta[tag] = float(meta[tag])

    meta['bands'] = {}
    mtag = 'Spectral_Band_Info'
    cdom =  xmldoc.getElementsByTagName(mtag)
    for cd in cdom:
        band = {}
        for tag in ['BAND_INDEX', 'BAND_DESCRIPTION', 'PHYSICAL_GAIN', 'PHYSICAL_BIAS', 'PHYSICAL_UNIT', 'ESUN']:
            tdom = cd.getElementsByTagName(tag)
            if len(tdom) > 0:
                if tdom[0].firstChild is None: continue
                band[tag] = tdom[0].firstChild.nodeValue
                if tag in ['PHYSICAL_GAIN', 'PHYSICAL_BIAS', 'ESUN']: band[tag] = float(band[tag])
                if tag in ['BAND_INDEX']:
                    band[tag] = int(band[tag])
                    band['INDEX'] = band[tag]-1

        if band['BAND_DESCRIPTION'] in ['BLUE', 'GREEN', 'RED']:
            band['BAND_DESCRIPTION'] = band['BAND_DESCRIPTION'].capitalize()
        meta['bands'][band['BAND_DESCRIPTION']] = band

    return(meta)
