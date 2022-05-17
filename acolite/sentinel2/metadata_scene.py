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
##                2021-10-13 (QV) adapted for new processing baseline which includes TOA offsets
##                2022-05-17 (QV) adapted for processing older N0201 files

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
    banddata['BandNames'] = {}
    banddata['Resolution'] = {}
    banddata['Wavelength'] = {}
    banddata['RSR'] = {}

    tdom = xmldoc.getElementsByTagName('Spectral_Information')
    for t in tdom:
        bandi = t.getAttribute('bandId')
        band = t.getAttribute('physicalBand')
        banddata['BandNames'][bandi] = band
        banddata['Resolution'][band] = t.getElementsByTagName('RESOLUTION')[0].firstChild.nodeValue
        banddata['Wavelength'][band] = {tag:float(t.getElementsByTagName(tag)[0].firstChild.nodeValue) for tag in ['CENTRAL','MIN','MAX']}
        tag = t.getElementsByTagName('Spectral_Response')
        if len(tag) > 0:
            step = float(tag[0].getElementsByTagName('STEP')[0].firstChild.nodeValue)
            rsr = [float(rs) for rs in tag[0].getElementsByTagName('VALUES')[0].firstChild.nodeValue.split(' ')]
            wave = np.linspace(banddata['Wavelength'][band]['MIN'],banddata['Wavelength'][band]['MAX'], int((banddata['Wavelength'][band]['MAX']-banddata['Wavelength'][band]['MIN'])/step)+1)
        banddata['RSR'][band] = {'response':rsr, 'wave':wave}

    ## workaround for N0201 data
    if len(banddata['BandNames']) == 0:
        tdom = xmldoc.getElementsByTagName('BAND_NAME')
        for ti, t in enumerate(tdom):
            bandi = '{}'.format(ti)
            band = t.firstChild.nodeValue
            banddata['BandNames'][bandi] = band
    ## if still empty
    if len(banddata['BandNames']) == 0:
        banddata['BandNames'] = {'0': 'B1', '1': 'B2', '2': 'B3', '3': 'B4', '4': 'B5', '5': 'B6', '6': 'B7',
                                 '7': 'B8', '8': 'B8A', '9': 'B9', '10': 'B10', '11': 'B11', '12': 'B12'}

    tdom = xmldoc.getElementsByTagName('SOLAR_IRRADIANCE')
    for t in tdom:
        if 'F0' not in banddata: banddata['F0'] = {}
        bandi = t.getAttribute('bandId')
        band = banddata['BandNames'][bandi]
        banddata['F0'][band] = float(t.firstChild.nodeValue) # 'unit':t.getAttribute('unit')

    tdom = xmldoc.getElementsByTagName('PHYSICAL_GAINS')
    for t in tdom:
        if 'PHYSICAL_GAINS' not in banddata: banddata['PHYSICAL_GAINS'] = {}
        bandi = t.getAttribute('bandId')
        band = banddata['BandNames'][bandi]
        banddata['PHYSICAL_GAINS'][band] = float(t.firstChild.nodeValue)

    tdom = xmldoc.getElementsByTagName('RADIO_ADD_OFFSET')
    for t in tdom:
        if 'RADIO_ADD_OFFSET' not in banddata: banddata['RADIO_ADD_OFFSET'] = {}
        bandi = t.getAttribute('band_id')
        band = banddata['BandNames'][bandi]
        banddata['RADIO_ADD_OFFSET'][band] = float(t.firstChild.nodeValue)

    return(metadata,banddata)
