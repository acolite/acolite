## def metadata_parse
## parses XML metadata from Pléiades images
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-06-29
## modifications: 2016-12-05 (QV) added check for wavelengths in nm and F0 == 999
##                2017-01-16 (QV) added read of extent vertices
##                2017-01-18 (QV) added pan option
##                2017-01-23 (QV) added band order information
##                2017-05-03 (QV) added reflectance gains if provided
##                2017-11-13 (QV) added computation of correct viewing azimuth angle
##                2017-11-21 (QV) added defaults for model selection and scene center geometry
##                2019-02-25 (QV) changed import name
##                2020-10-18 (QV) added tile info
##                2020-10-18 (QV) renamed from parse_metadata, prep for acg


def metadata_parse(metafile, pan=False):
    import os, sys, fnmatch, dateutil.parser
    from xml.dom import minidom
    import numpy as np
    dtor = np.pi/180.

    xmldoc = minidom.parse(metafile)

    metadata = {}
    metadata['satellite']='Pléiades'

    ## get image information
    metadata_tags = ['NROWS','NCOLS','NBANDS', 'RESAMPLING_SPACING', 'MISSION','MISSION_INDEX',
                     'INSTRUMENT','INSTRUMENT_INDEX','IMAGING_DATE', 'IMAGING_TIME', 'BAND_MODE',
                     'RED_CHANNEL', 'GREEN_CHANNEL', 'BLUE_CHANNEL', 'ALPHA_CHANNEL', 'EXTENT_TYPE']
    for tag in metadata_tags:
        node = xmldoc.getElementsByTagName(tag)
        if len(node) > 0: metadata[tag] = node[0].firstChild.nodeValue
    metadata['sensor'] = '{}{}'.format(metadata['MISSION'],metadata['MISSION_INDEX'])
    if 'SPOT' in metadata['sensor']: metadata['satellite'] = 'SPOT'

    ## get acquisition time
    metadata["isotime"] = '{}T{}'.format(metadata["IMAGING_DATE"],metadata["IMAGING_TIME"])

    ###
    node = xmldoc.getElementsByTagName("RADIOMETRIC_PROCESSING")
    metadata['RADIOMETRIC_PROCESSING'] = node[0].firstChild.nodeValue if len(node) > 0 else 'RADIANCE'

    ## get 'special' values
    for t in xmldoc.getElementsByTagName('Special_Value'):
        name = (t.getElementsByTagName("SPECIAL_VALUE_TEXT")[0].firstChild.nodeValue)
        value = (t.getElementsByTagName("SPECIAL_VALUE_COUNT")[0].firstChild.nodeValue)
        metadata[name] = value

    ## get view and sun geometry
    geometric_tags = ['LOCATION_TYPE','TIME', 'SUN_AZIMUTH','SUN_ELEVATION',
                      'AZIMUTH_ANGLE','VIEWING_ANGLE_ACROSS_TRACK','VIEWING_ANGLE_ALONG_TRACK','VIEWING_ANGLE',
                      'INCIDENCE_ANGLE_ALONG_TRACK','INCIDENCE_ANGLE_ACROSS_TRACK','INCIDENCE_ANGLE']
    geometric_values = []
    for t in xmldoc.getElementsByTagName('Located_Geometric_Values') :
        geom={}
        for tag in geometric_tags:
            node = t.getElementsByTagName(tag)
            if len(node) > 0:
                if tag in ['LOCATION_TYPE','TIME']: geom[tag]=node[0].firstChild.nodeValue
                else: geom[tag]=float(node[0].firstChild.nodeValue)

        ## compute viewing azimuth
        orientation_angle = geom['AZIMUTH_ANGLE']
        ang_across = geom['INCIDENCE_ANGLE_ACROSS_TRACK']
        ang_along = geom['INCIDENCE_ANGLE_ALONG_TRACK']
        geom['VIEWING_AZIMUTH'] = np.mod(orientation_angle - np.arctan2(np.tan(ang_across*dtor), np.tan(ang_along*dtor))/dtor,360)

        geometric_values.append(geom)

    metadata['GEOMETRY'] = geometric_values

    ## get band info such as F0 and calibration
    band_info = {}
    bands = ['B0','B1','B2','B3']
    default_F0 = [1915., 1830., 1594.,1060.]

    if pan is True:
        bands = ['P']
        default_F0 = [1548.]

    for band in bands:
        band_data = {}
        for t in xmldoc.getElementsByTagName('BAND_ID') :
            if (t.firstChild.nodeValue) == band:
                parent = t.parentNode.nodeName
                if parent == 'Band_Solar_Irradiance':
                    unit = t.parentNode.getElementsByTagName('MEASURE_UNIT')[0].firstChild.nodeValue
                    F0 = float(t.parentNode.getElementsByTagName('VALUE')[0].firstChild.nodeValue)
                    if F0 == 999.:
                         idx = [i for i,j in enumerate(bands) if j == band]
                         F0 = default_F0[idx[0]]
                    band_data['F0'] = F0
                if parent == 'Band_Radiance':
                    unit = t.parentNode.getElementsByTagName('MEASURE_UNIT')[0].firstChild.nodeValue
                    band_data['radiance_gain'] = float(t.parentNode.getElementsByTagName('GAIN')[0].firstChild.nodeValue)
                    band_data['radiance_bias'] = float(t.parentNode.getElementsByTagName('BIAS')[0].firstChild.nodeValue)

                if parent == 'Band_Reflectance':
                    band_data['reflectance_gain'] = float(t.parentNode.getElementsByTagName('GAIN')[0].firstChild.nodeValue)
                    band_data['reflectance_bias'] = float(t.parentNode.getElementsByTagName('BIAS')[0].firstChild.nodeValue)

                if parent == 'Band_Digital_Number':
                    band_data['radiance_to_dn_gain'] = float(t.parentNode.getElementsByTagName('GAIN')[0].firstChild.nodeValue)
                    band_data['radiance_to_dn_bias'] = float(t.parentNode.getElementsByTagName('BIAS')[0].firstChild.nodeValue)

                if parent == 'Band_Spectral_Range':
                    unit = t.parentNode.getElementsByTagName('MEASURE_UNIT')[0].firstChild.nodeValue
                    wave_start = float(t.parentNode.getElementsByTagName('MIN')[0].firstChild.nodeValue)
                    wave_end = float(t.parentNode.getElementsByTagName('MAX')[0].firstChild.nodeValue)
                    if wave_start > 100.: wave_start/=1000.
                    if wave_end > 100.: wave_end/=1000.
                    band_data['wave_start'] = wave_start
                    band_data['wave_end'] = wave_end
                    band_data['wave'] = (band_data['wave_end']+band_data['wave_start'])/2.
                    band_data['wave_name'] = str(round(int(band_data['wave']*1000.),2))
                if pan is False:
                    if band == metadata['RED_CHANNEL']: band_data['band_index']=0
                    if band == metadata['GREEN_CHANNEL']: band_data['band_index']=1
                    if band == metadata['BLUE_CHANNEL']: band_data['band_index']=2
                    if band == metadata['ALPHA_CHANNEL']: band_data['band_index']=3
        band_info[band]=band_data
    metadata['BAND_INFO'] = band_info

    ## read vertices of spatial extent
    vertex_tags = ['LON','LAT','X','Y','COL','ROW']
    vertex_values=[]
    for t in xmldoc.getElementsByTagName('Vertex'):
        vertex = {}
        for tag in vertex_tags:
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    vertex[tag]=float(node[0].firstChild.nodeValue)
        vertex_values.append(vertex)

    vertices = {}
    for t in xmldoc.getElementsByTagName('Center'):
        vertex = {}
        for tag in vertex_tags:
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    vertex[tag]=float(node[0].firstChild.nodeValue)
        vertices['C']=vertex

    ncols = float(metadata['NCOLS'])
    nrows = float(metadata['NROWS'])
    for i,v in enumerate(vertex_values):
        col = float(v['COL'])
        row = float(v['ROW'])
        if (col == 1) & (row == 1):
            vertices['UL'] = v
            continue
        if (col == ncols) & (row == 1):
            vertices['UR'] = v
            continue
        if (col == ncols) & (row == nrows):
            vertices['LR'] = v
            continue
        if (col == 1) & (row == nrows):
            vertices['LL'] = v
            continue
        vertices['V{}'.format(i)]=v

    metadata['VERTICES']=vertices

    ## get tile information
    try:
        metadata['ntiles'] = int(xmldoc.getElementsByTagName("NTILES")[0].firstChild.nodeValue)
        metadata['ntiles_R'] = int(xmldoc.getElementsByTagName("NTILES_COUNT")[0]._attrs['ntiles_R'].nodeValue)
        metadata['ntiles_C'] = int(xmldoc.getElementsByTagName("NTILES_COUNT")[0]._attrs['ntiles_C'].nodeValue)
        metadata['tiles_nrows'] = int(xmldoc.getElementsByTagName("NTILES_SIZE")[0]._attrs['nrows'].nodeValue)
        metadata['tiles_ncols'] = int(xmldoc.getElementsByTagName("NTILES_SIZE")[0]._attrs['ncols'].nodeValue)
    except:
        metadata['ntiles'] = 1
        metadata['ntiles_R'] = 1
        metadata['ntiles_C'] = 1
        metadata['tiles_nrows'] = int(metadata['NROWS'])
        metadata['tiles_ncols'] = int(metadata['NROWS'])

    ## set some defaults
    metadata['sza'] = 90. - metadata['GEOMETRY'][1]['SUN_ELEVATION']
    metadata['vza'] = metadata['GEOMETRY'][1]['VIEWING_ANGLE']
    metadata['raa'] = abs(metadata['GEOMETRY'][1]['SUN_AZIMUTH'] - metadata['GEOMETRY'][1]['VIEWING_AZIMUTH'])
    while metadata['raa'] > 180: metadata['raa'] = abs(metadata['raa']-360)

    return(metadata)
