## def metadata_parse
## parses XML metadata from PlanetScope bundle images
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-12
## modifications: 2018-03-15 (QV) moved band stuff to separate keys
##                2018-03-27 (QV) added 0d = 0c
##                2019-06-12 (QV) added 11 = 0f
##                2021-02-24 (QV) reformatted and simplified for acg
##                2022-02-21 (QV) added Skysat

def metadata_parse(metafile):
    import os, json
    from xml.dom import minidom
    import dateutil.parser

    fn, ext = os.path.splitext(metafile)
    if ext == '.json': ## json metadata
        with open(metafile) as f:
            metadata = json.load(f)

        try:
            if metadata['properties']['provider'] == 'skysat': metadata['platform'] = 'Skysat'
        except:
            print('Platform not recognised {}'.format(metafile))

        metadata['satellite_id'] = metadata['properties']['satellite_id']
        metadata['isotime'] = metadata['properties']['acquired']
        metadata['sensor'] = metadata['properties']['satellite_id'].replace('SSC', 'Skysat')
        metadata['satellite_sensor'] = metadata['properties']['satellite_id'].replace('SSC', 'Skysat')
        if metadata['satellite_sensor'] in ['Skysat14', 'Skysat15', 'Skysat16', 'Skysat17', 'Skysat18', 'Skysat19']:
            metadata['sensor'] = 'Skysat14-SkySat19'

        metadata['sza'] = 90-metadata['properties']['sun_elevation']
        metadata['saa'] = metadata['properties']['sun_azimuth']
        metadata['vza'] = metadata['properties']['view_angle']
        metadata['vaa'] = metadata['properties']['satellite_azimuth']
        metadata['raa'] = abs(metadata['saa'] - metadata['vaa'])
        while metadata['raa'] > 180: metadata['raa'] = abs(metadata['raa']-360)

        metadata['Blue-band_idx'] = 1
        metadata['Green-band_idx'] = 2
        metadata['Red-band_idx'] = 3
        metadata['NIR-band_idx'] = 4

    elif ext == '.xml': ## xml metadata
        xmldoc = minidom.parse(metafile)
        metadata = {}

        ## get platform info
        main_tag = 'eop:Platform'
        tags = ["eop:shortName", "eop:serialIdentifier","eop:orbitType"]
        tags_out = ["platform", 'platform_id', 'orbit']
        for t in xmldoc.getElementsByTagName('eop:Platform'):
            for i,tag in enumerate(tags):
                    node = t.getElementsByTagName(tag)
                    if len(node) > 0:
                        metadata[tags_out[i]] = node[0].firstChild.nodeValue

        if "RapidEye" in metadata['platform_id']:
            metadata['satellite_id'] = metadata['platform']
            metadata['satellite_sensor'] = metadata['platform_id']
            metadata['sensor'] = "RapidEye"
            metadata['sensor_family']='RapidEye'
            metadata['satellite_prefix'] = 're'

        if 'PlanetScope' in metadata['platform']:
            metadata['satellite_id'] = metadata['platform_id'][0:2]
            if metadata['satellite_id'] == '10': metadata['satellite_id'] = "0f"
            if metadata['satellite_id'] == '11': metadata['satellite_id'] = "0f"
            if metadata['satellite_id'] == '0d': metadata['satellite_id'] = "0c"

            ## assume 23/24XX doves use the same RSR as the 22XX
            if metadata['satellite_id'] >= '22': metadata['satellite_id'] = "22"
            metadata['sensor_family']='PlanetScope'
            metadata['satellite_prefix'] = 'ps'

        ## get band data
        main_tag='{}:bandSpecificMetadata'.format(metadata['satellite_prefix'])
        bands0 = {}
        tags = ["{}:bandNumber".format(metadata['satellite_prefix']),
                '{}:radiometricScaleFactor'.format(metadata['satellite_prefix']),
                '{}:reflectanceCoefficient'.format(metadata['satellite_prefix'])]
        tags_out = ["band_idx",'to_radiance', 'to_reflectance']
        for t in xmldoc.getElementsByTagName(main_tag):
            band = {}
            for i,tag in enumerate(tags):
                    node = t.getElementsByTagName(tag)
                    if len(node) > 0:
                        val = float(node[0].firstChild.nodeValue)
                        band[tags_out[i]] = val
            bands0[band['band_idx']] = band

        bnames = {}
        if metadata['satellite_prefix'] == 'ps':
            ## default 4B data
            if (len(bands0) == 4):
                bnames = {1:'Blue',2:'Green',3:'Red',4:'NIR'}

            if (metadata['satellite_id'] == '22') & (len(bands0) == 4):
                bnames = {1:'Blue',2:'Green',3:'Red',4:'NIR'}
            if (metadata['satellite_id'] == '22') & (len(bands0) == 5):
                bnames = {1:'Blue',2:'Green_ii',3:'Red',4:'RedEdge', 5:'NIR'}
                metadata['satellite_id'] = 'SD5'
            if (metadata['satellite_id'] == '22') & (len(bands0) == 8):
                bnames = { 1:'Coastal-Blue',  2:'Blue', 3:'Green_i', 4:'Green_ii',
                           5:'Yellow', 6:'Red', 7:'Red-edge', 8:'NIR'}
                metadata['satellite_id'] = 'SD8'

            metadata['satellite_sensor'] = '{}_{}'.format(metadata['platform'], metadata['platform_id'])
            metadata['sensor'] = '{}_{}'.format(metadata['platform'], metadata['satellite_id'])

        if metadata['satellite_prefix'] == 're':
            bnames = {1:'Blue',2:'Green',3:'Red',4:'RedEdge', 5:'NIR'}

        if len(bnames) == 0:
            print('Error determining bands from metadata')

        bands = {bnames[bands0[b]['band_idx']]:bands0[b] for b in bands0}
        for band in bands:
            for key in bands[band]:
                bk='{}-{}'.format(band, key)
                metadata[bk] = bands[band][key]

        ## get acquisition info
        main_tag = 'eop:acquisitionParameters'
        tags = ["eop:orbitDirection", "eop:incidenceAngle",
                "opt:illuminationAzimuthAngle", "opt:illuminationElevationAngle",
                "{}:azimuthAngle".format(metadata['satellite_prefix']),
                "{}:spaceCraftViewAngle".format(metadata['satellite_prefix']),
                "{}:acquisitionDateTime".format(metadata['satellite_prefix'])]
        tags_out = ["orbit",'ViewingIncidence',
                    'SunAzimuth', 'SunElevation',
                    'ViewingAzimuth', 'ViewZenith', 'isotime']
        for t in xmldoc.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                    node = t.getElementsByTagName(tag)
                    if len(node) > 0:
                        if tag in ["eop:orbitDirection","{}:acquisitionDateTime".format(metadata['satellite_prefix'])]:
                            val=node[0].firstChild.nodeValue
                        else:
                            val=float(node[0].firstChild.nodeValue)
                        metadata[tags_out[i]] = val

        ## geometry info
        metadata['sza'] = 90.-metadata['SunElevation']
        metadata['saa'] = metadata['SunAzimuth']
        metadata['vza'] = abs(metadata['ViewZenith'])
        metadata['vaa'] = metadata['ViewingAzimuth']
        metadata['raa'] = abs(metadata['saa'] - metadata['vaa'])
        while metadata['raa'] > 180: metadata['raa'] = abs(metadata['raa']-360)

        ## get product info
        main_tag = '{}:spatialReferenceSystem'.format(metadata['satellite_prefix'])
        tags = ["{}:epsgCode".format(metadata['satellite_prefix']),
                "{}:geodeticDatum".format(metadata['satellite_prefix']),
                "{}:projection".format(metadata['satellite_prefix']),
                "{}:projectionZone".format(metadata['satellite_prefix'])]
        tags_out = ["epsg",'datum', 'projection', 'zone']
        for t in xmldoc.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                    node = t.getElementsByTagName(tag)
                    if len(node) > 0:
                        val=node[0].firstChild.nodeValue
                        metadata[tags_out[i]] = val

        ## get resolution info
        main_tag = '{}:ProductInformation'.format(metadata['satellite_prefix'])
        tags = ["{}:numRows".format(metadata['satellite_prefix']),
                "{}:numColumns".format(metadata['satellite_prefix']),
                "{}:numBands".format(metadata['satellite_prefix']),
                "{}:rowGsd".format(metadata['satellite_prefix']),
                "{}:columnGsd".format(metadata['satellite_prefix'])]
        tags_out = ["nrow","ncol","nband","resolution_row", "resolution_col"]
        for t in xmldoc.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                    node = t.getElementsByTagName(tag)
                    if len(node) > 0:
                        val=node[0].firstChild.nodeValue
                        metadata[tags_out[i]] = val
        metadata['resolution'] = (float(metadata['resolution_row']),float(metadata['resolution_col']))
        metadata['dims'] = (int(metadata['ncol']),int(metadata['nrow']))

        ## get bounding box
        main_tag = '{}:geographicLocation'.format(metadata['satellite_prefix'])
        tags = ["{}:topLeft".format(metadata['satellite_prefix']),
                "{}:topRight".format(metadata['satellite_prefix']),
                "{}:bottomRight".format(metadata['satellite_prefix']),
                "{}:bottomLeft".format(metadata['satellite_prefix'])]
        tags_out = ["UL",'UR', 'LR', 'LL']
        for t in xmldoc.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                    node = t.getElementsByTagName(tag)
                    for j,tag2 in enumerate(['{}:latitude'.format(metadata['satellite_prefix']),'{}:longitude'.format(metadata['satellite_prefix'])]):
                        node2 = node[0].getElementsByTagName(tag2)
                        if len(node2) > 0:
                            val=node2[0].firstChild.nodeValue
                            tout = '{}_{}'.format(tags_out[i], tag2.split(':')[1].upper())
                            metadata[tout] = float(val)
    return(metadata)
