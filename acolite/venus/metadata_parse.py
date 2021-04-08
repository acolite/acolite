## def metadata_parse
## parses XML metadata from VENUS images
## written by Quinten Vanhellemont, RBINS
## 2021-04-08
## modifications:

def metadata_parse(bundle):
    import glob
    from xml.dom import minidom
    import numpy as np

    ## find metadata file
    metafile = glob.glob('{}/*.xml'.format(bundle))
    if len(metafile) == 1:
        metafile = metafile[0]
    else: return()

    print(metafile)
    ## open xml
    xmldoc = minidom.parse(metafile)

    metadata = {}

    ## get image information
    metadata_tags = ['PLATFORM', 'SPECTRAL_CONTENT', 'REFLECTANCE_QUANTIFICATION_VALUE']
    for tag in metadata_tags:
        node = xmldoc.getElementsByTagName(tag)
        if len(node) > 0: metadata[tag] = node[0].firstChild.nodeValue
    metadata['REFLECTANCE_QUANTIFICATION_VALUE'] = float(metadata['REFLECTANCE_QUANTIFICATION_VALUE'])

    characteristics = xmldoc.getElementsByTagName('Product_Characteristics')[0]
    metadata_tags = ['ACQUISITION_DATE']
    for tag in metadata_tags:
        node = characteristics.getElementsByTagName(tag)
        if len(node) > 0: metadata[tag] = node[0].firstChild.nodeValue

    ## image properties
    im = xmldoc.getElementsByTagName('Image_Properties')[0]
    for k in ['NATURE', 'FORMAT', 'ENCODING', 'ENDIANNESS']:
        #print(k, im.getElementsByTagName(k)[0].firstChild.nodeValue)
        if k == 'NATURE':
            metadata['image_type'] = im.getElementsByTagName(k)[0].firstChild.nodeValue

    ## get geolocation information
    geolocation = {'points': {}}
    gg = xmldoc.getElementsByTagName('Global_Geopositioning')[0]
    for n in gg.getElementsByTagName('Point'):
        name = n.getAttribute('name')
        geolocation['points'][name] = {}
        for k in ['LAT', 'LON', 'X', 'Y']:
            for p in n.getElementsByTagName(k):
                geolocation['points'][name][k] = float(p.firstChild.nodeValue)
    for k in ['GEO_TABLES', 'HORIZONTAL_CS_TYPE', 'HORIZONTAL_CS_NAME', 'HORIZONTAL_CS_CODE']:
        cur = xmldoc.getElementsByTagName(k)[0]
        geolocation[k] = cur.firstChild.nodeValue

    metadata['geolocation'] = geolocation
    metadata['xrange'] = metadata['geolocation']['points']['upperLeft']['X'], metadata['geolocation']['points']['lowerRight']['X']
    metadata['yrange'] = metadata['geolocation']['points']['upperLeft']['Y'], metadata['geolocation']['points']['lowerRight']['Y']

    ## scene dimensions
    metadata_tags = ['ULX', 'ULY', 'XDIM', 'YDIM', 'NROWS', 'NCOLS']
    for tag in metadata_tags:
        node = xmldoc.getElementsByTagName(tag)
        if len(node) > 0: metadata[tag] = float(node[0].firstChild.nodeValue)

    #metadata['xrange'] = metadata['ULX'] + metadata['XDIM'] * metadata['NCOLS']
    #metadata['xrange'] = metadata['ULY'] + metadata['YDIM'] * metadata['NROWS']

    ## get bands information
    bands = {}
    masks = {}
    ## bands
    for n in xmldoc.getElementsByTagName('IMAGE_FILE'):
        band = n.getAttribute('band_id')
        image_file = n.firstChild.nodeValue
        image_path = '{}/{}'.format(bundle, image_file)
        bands[band] = {'file':image_file, 'path':image_path}
    ## masks
    for n in xmldoc.getElementsByTagName('MASK_FILE'):
        band = n.getAttribute('band_id')
        mask_file = n.firstChild.nodeValue
        mask_path = '{}/{}'.format(bundle, mask_file)
        mask_type = None
        for k in ['_PIX_', '_SAT_', '_CLD_', '_USI_']:
            if k in mask_file: mask_type = k.strip('_')
        if band in bands:
            bands[band][mask_type] = {'file':mask_file, 'path':mask_path}
        else:
            masks[mask_type] = {'file':mask_file, 'path':mask_path}
    metadata['bands'] = bands
    metadata['masks'] = masks

    ## get additional data info
    data_list = []
    tmp = xmldoc.getElementsByTagName('Data_List')[0].getElementsByTagName('Data')
    for m in tmp:
        for k in ['Data_Properties']:
            name = m.getElementsByTagName(k)[0].getElementsByTagName('NATURE')[0].firstChild.nodeValue
        for n in m.getElementsByTagName('DATA_FILE'):
            properties = {'name':name}
            fname = n.firstChild.nodeValue
            properties['file'] = fname
            properties['path'] = '{}/{}'.format(bundle, fname)
            for a in ['altitude', 'axis', 'band_number', 'group_id']:
                att = n.getAttribute(a)
                if len(att) == 0: continue
                properties[a] = att
            data_list.append(properties)
    metadata['data'] = data_list

    ## sun/view geometry average
    tmp = xmldoc.getElementsByTagName('Mean_Value_List')[0]
    sun = tmp.getElementsByTagName('Sun_Angles')[0]

    sza = float(sun.getElementsByTagName('ZENITH_ANGLE')[0].firstChild.nodeValue)
    saa = float(sun.getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.nodeValue)

    tmp = xmldoc.getElementsByTagName('Mean_Viewing_Incidence_Angle_List')[0]
    view = tmp.getElementsByTagName('Mean_Viewing_Incidence_Angle')
    view_angles = {}

    for m in view:
        det = m.getAttribute('detector_id')
        vza = m.getElementsByTagName('ZENITH_ANGLE')[0].firstChild.nodeValue
        vaa = m.getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.nodeValue
        view_angles[det] = {'vza':float(vza), 'vaa':float(vaa)}

    vza = np.nanmean([view_angles[d]['vza'] for d in view_angles])
    vaa = np.nanmean([view_angles[d]['vaa'] for d in view_angles])
    metadata['vza'] = vza
    metadata['vaa'] = vaa
    metadata['sza'] = sza
    metadata['saa'] = saa
    metadata['raa'] = saa-vaa
    if metadata['raa'] > 180:
        metadata['raa'] -= 180
    if metadata['raa'] < 0:
        metadata['raa'] = np.abs(metadata['raa'])

    ## set some common variables
    if metadata['PLATFORM'] == 'VENUS':
        metadata['sensor'] = 'VENµS_VSSC'
    metadata['isodate'] = metadata['ACQUISITION_DATE']
    metadata['pixel_size'] = metadata['XDIM'], metadata['YDIM']
    metadata['scene_dims'] = metadata['NCOLS'], metadata['NROWS']

    ## projection info
    cs_code = metadata['geolocation']['HORIZONTAL_CS_CODE']
    epsg = int(cs_code)
    datum = 'WGS84'
    if 32600 < epsg <= 32660:
        zone = epsg - 32600
        proj4_list = ['+proj=utm',
                      '+zone={}'.format(zone),
                      '+datum={}'.format(datum),
                      '+units=m',
                      '+no_defs ']
    if 32700 < epsg <= 32760:
        zone = epsg - 32700
        proj4_list = ['+proj=utm',
                      '+zone={}'.format(zone),
                      '+south',
                      '+datum={}'.format(datum),
                      '+units=m',
                      '+no_defs ']
    metadata['zone'] = zone
    metadata['proj4_string'] = ' '.join(proj4_list)

    return(metadata)
