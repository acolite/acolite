## def metadata_parse
## parses XML metadata from DESIS bundle images
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-27
## edited by Dirk Aurin, NASA Goddard Space Flight Center
## modifications: QV 2019-05-06
##                2020-03-16 (QV) fixed tile metadata reads
##                2020-07-26 (QV) added WV3
##                2021-02-24 (QV) reformatted and simplified for acg
##                2021-06-28 (DA) adapt from WV2/WV3 to DESIS

def metadata(metafile):
    import os, sys, fnmatch, dateutil.parser
    from xml.dom import minidom
    from osgeo import gdal

    if not os.path.isfile(metafile):
        print('Metadata file not found.')
        sys.exit()

    try:
        xmldoc = minidom.parse(metafile)
    except:
        print('Error opening metadata file.')
        sys.exit()

    metadata = {}
    ## get image information    
    metadata_tags = ['backgroundValue','version','l1aVersion','mission', 'satelliteID', 'sensor', 
                        'tileID', 'startTime', 'endTime','productType', 'numberOfBands', 
                        'orbitDirection', 'sceneAzimuthAngle', 'sceneIncidenceAngle', 
                        'sunAzimuthAngle', 'sunZenithAngle', 'numberOfBands']

    for tag in metadata_tags:
        node = xmldoc.getElementsByTagName(tag)
        if len(node) > 0: metadata[tag] = node[0].firstChild.nodeValue

    ## Get spatial coverage (only necessary for L1B)
    ## Read in point info
    point_tags = ['latitude','longitude']
    point_values={}
        
    for i, t in enumerate(xmldoc.getElementsByTagName('point')):
        node = t.getElementsByTagName('frame')
        if len(node) > 0:
            point_name = node[0].firstChild.nodeValue
            point_data = {'name': point_name}
            point_data['index'] = i+1    # Not zero referenced??

            for tag in point_tags:
                    node = t.getElementsByTagName(tag)
                    if len(node) > 0:
                        # band_data[tag]=float(node[0].firstChild.nodeValue)
                        value = node[0].firstChild.nodeValue
                        valList = value.split(',')
                        if len(valList)==1:
                            point_data[tag]=float(value)
                        else:    
                            point_data[tag]=[float(x) for x in  valList]
                            
        if len(point_data)>2: ## keep point only if in metadata
            point_values[point_name] = point_data


    # Raster dimensions are not included in DESIS metadata xml, so read in the TIF
    fpList = os.path.split(metafile)
    file = '{}/{}'.format(fpList[0],fpList[1].replace('METADATA.xml','SPECTRAL_IMAGE.tif'))
    ## check if the files were named .TIF instead of .tif
    if not os.path.exists(file): file = file.replace('.tif', '.TIF')
    
    try:
        ## open file
        ds = gdal.Open(file)        
        dimx, dimy = ds.RasterXSize, ds.RasterYSize
        ds = None
    except:
        print('Could not determine projection from {}'.format(file))
        pass

    if 'mission' in metadata:
        if metadata['mission'] == 'DESIS':
            metadata['satellite'] = metadata['satelliteID']
            # metadata['sensor'] = f"{metadata['mission']}_{metadata['sensor']}"
            metadata['isotime']=metadata['startTime'] #metadata["EARLIESTACQTIME"]
            
            ## geometry info
            metadata['sza'] = float(metadata['sunZenithAngle'])
            metadata['saa'] = float(metadata['sunAzimuthAngle'])
            metadata['vza'] = float(metadata['sceneIncidenceAngle'])
            metadata['vaa'] = float(metadata['sceneAzimuthAngle'])
            metadata['raa'] = abs(metadata['saa'] - metadata['vaa'])
            while metadata['raa'] > 180: metadata['raa'] = abs(metadata['raa']-360)
            metadata['NUMROWS'] = dimy
            metadata['NUMCOLUMNS'] = dimx    

            ## Interpretation of points depends on orbit direction
            #   "The bounding polygon is a list of points starting with the most starboard pixel 
            #   of the first scanline and followed by a counter-clockwise sequence of four points,
            #   whereas the last point is identical with the first point."
            # Assuming "Left" and "Uppoer", etc. refer to how the image was collected, this should 
            # translate to point_1: LR, 2: UR, 3: UL, 4: LL 
            #   However, it is not clear that this is how the 1024x1024 image in L1B is oriented. They 
            #   ARE tilted from north, but are they ever upside down (i.e. descending)?
            # This should not be necessary using L1C orthrectified and projected images. These would 
            # only be used in the case that Gdal is unable to read the projection info in the GeoTiff, which is 
            # hit-and-miss at L1B, apparently regardless of processing version. 
            metadata['LLLON'] = point_values['point_4']['longitude']       
            metadata['LLLAT'] = point_values['point_4']['latitude'] 
            metadata['ULLON'] = point_values['point_3']['longitude'] 
            metadata['ULLAT'] = point_values['point_3']['latitude'] 
            metadata['URLON'] = point_values['point_2']['longitude'] 
            metadata['URLAT'] = point_values['point_2']['latitude'] 
            metadata['LRLON'] = point_values['point_1']['longitude'] 
            metadata['LRLAT'] = point_values['point_1']['latitude'] 

    ## Read in band info
    band_tags = ['wavelengthCenterOfBand', 'wavelengthWidthOfBand', 'response',
                'wavelengths', 'gainOfBand', 'offsetOfBand']
    band_values={}
        
    for i, t in enumerate(xmldoc.getElementsByTagName('band')):
        node = t.getElementsByTagName('bandNumber')
        if len(node) > 0:
            band_name = node[0].firstChild.nodeValue
            band_data = {'name': band_name}
            band_data['index'] = i+1    # Not zero referenced??

            for tag in band_tags:
                    node = t.getElementsByTagName(tag)
                    if len(node) > 0:
                        # band_data[tag]=float(node[0].firstChild.nodeValue)
                        value = node[0].firstChild.nodeValue
                        valList = value.split(',')
                        if len(valList)==1:
                            band_data[tag]=float(value)
                        else:    
                            band_data[tag]=[float(x) for x in  valList]
                            
        if len(band_data)>2: ## keep band only if in metadata
            band_values[band_name] = band_data

    metadata['BAND_INFO'] = band_values


    # ## get tile information
    # tile_tags = ["FILENAME",
    #             "ULCOLOFFSET","ULROWOFFSET","URCOLOFFSET","URROWOFFSET",
    #             "LRCOLOFFSET","LRROWOFFSET","LLCOLOFFSET","LLROWOFFSET",
    #             "ULLON","ULLAT","URLON","URLAT",
    #             "LRLON","LRLAT","LLLON","LLLAT",
    #             "ULX","ULY","URX","URY",
    #             "LRX","LRY","LLX","LLY"]
    # tile_values=[]
    # for t in xmldoc.getElementsByTagName('TILE'):
    #     tile = {}
    #     for tag in tile_tags:
    #             node = t.getElementsByTagName(tag)
    #             if len(node) > 0:
    #                 if tag == "FILENAME": val=node[0].firstChild.nodeValue
    #                 else: val=float(node[0].firstChild.nodeValue)
    #                 tile[tag]=val
    #     tile_values.append(tile)
    # metadata['TILE_INFO'] = tile_values

    return(metadata)
