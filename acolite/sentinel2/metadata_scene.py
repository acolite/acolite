## def metadata_scene
## imports S2 scene metadata
## written by Quinten Vanhellemont, RBINS
## 2017-04-18
## modifications: 2017-04-27 (QV) added bandnames fallback for older metadata
##                2017-06-06 (QV) added waves and rgb bands info
##                2017-11-22 (QV) added defaults for the model selection
##                2018-07-18 (QV) changed acolite import name
##                2018-10-01 (QV) removed obsolete bits
##                2021-02-11 (QV) adapted for acolite-gen, renamed from scene_meta


def metadata_scene(metafile):
    import dateutil.parser
    from xml.dom import minidom

    #from acolite.shared import distance_se
    import numpy as np

    try:
        xmldoc = minidom.parse(metafile)
    except:
        print('Error opening metadata file.')
        sys.exit()

    xml_main = xmldoc.firstChild

    metadata = {}

    tags = ['PRODUCT_START_TIME','PRODUCT_STOP_TIME','PRODUCT_URI','PROCESSING_LEVEL',
            'PRODUCT_TYPE', 'PROCESSING_BASELINE', 'GENERATION_TIME','SPACECRAFT_NAME',
            'DATATAKE_SENSING_START', 'SENSING_ORBIT_NUMBER', 'SENSING_ORBIT_DIRECTION',
            'PRODUCT_FORMAT', 'QUANTIFICATION_VALUE', 'U']

    for tag in tags:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) > 0:
            if tdom[0].firstChild is not None:
                metadata[tag] = tdom[0].firstChild.nodeValue

    ##
    tdom = xmldoc.getElementsByTagName('Special_Values')
    for t in tdom:
        fill = (t.getElementsByTagName('SPECIAL_VALUE_TEXT')[0].firstChild.nodeValue)
        fill_value = (t.getElementsByTagName('SPECIAL_VALUE_INDEX')[0].firstChild.nodeValue)
        metadata[fill] = fill_value

    ## get information for sensor bands
    banddata = {}
    banddata['F0'] = {}
    tdom = xmldoc.getElementsByTagName('SOLAR_IRRADIANCE')
    for t in tdom:
        band = t.getAttribute('bandId')
        banddata['F0'][band] = float(t.firstChild.nodeValue) # 'unit':t.getAttribute('unit')

    banddata['PHYSICAL_GAINS'] = {}
    tdom = xmldoc.getElementsByTagName('PHYSICAL_GAINS')
    for t in tdom:
        band = t.getAttribute('bandId')
        banddata['PHYSICAL_GAINS'][band] = float(t.firstChild.nodeValue)

    banddata['BandNames'] = {}
    banddata['Resolution'] = {}
    banddata['Wavelength'] = {}
    banddata['RSR'] = {}

    tdom = xmldoc.getElementsByTagName('Spectral_Information')
    for t in tdom:
        band = t.getAttribute('bandId')
        banddata['BandNames'][band] = t.getAttribute('physicalBand')
        banddata['Resolution'][band] = t.getElementsByTagName('RESOLUTION')[0].firstChild.nodeValue
        banddata['Wavelength'][band] = {tag:float(t.getElementsByTagName(tag)[0].firstChild.nodeValue) for tag in ['CENTRAL','MIN','MAX']}
        tag = t.getElementsByTagName('Spectral_Response')
        if len(tag) > 0:
            step = float(tag[0].getElementsByTagName('STEP')[0].firstChild.nodeValue)
            rsr = [float(rs) for rs in tag[0].getElementsByTagName('VALUES')[0].firstChild.nodeValue.split(' ')]
            wave = np.linspace(banddata['Wavelength'][band]['MIN'],banddata['Wavelength'][band]['MAX'], int((banddata['Wavelength'][band]['MAX']-banddata['Wavelength'][band]['MIN'])/step)+1)
        banddata['RSR'][band] = {'response':rsr, 'wave':wave}

    ## some interpretation
    #metadata['TIME'] = dateutil.parser.parse(metadata['PRODUCT_STOP_TIME'])
    #metadata["DOY"] = metadata["TIME"].strftime('%j')
    #metadata["SE_DISTANCE"] = distance_se(metadata['DOY'])

    #metadata["SATELLITE"] = metadata['SPACECRAFT_NAME'] # 'Sentinel-2A'

    return(metadata,banddata)
