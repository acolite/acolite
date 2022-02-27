## def metadata_parse
## parses XML metadata from WorldView2 bundle images
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-27
## modifications: QV 2019-05-06
##                2020-03-16 (QV) fixed tile metadata reads
##                2020-07-26 (QV) added WV3
##                2021-02-24 (QV) reformatted and simplified for acg
##                2022-02-24 (QV) added GeoEye1
##                2022-02-27 (QV) change band_index according to band tags present in metadata

def metadata_parse(metafile):
    import os, sys, fnmatch, dateutil.parser
    from xml.dom import minidom

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
    metadata_tags = ['SATID', 'FIRSTLINETIME', 'NUMROWS','NUMCOLUMNS','PRODUCTLEVEL'
                    "MININTRACKVIEWANGLE", "MAXINTRACKVIEWANGLE", "MEANINTRACKVIEWANGLE",
                    "MINCROSSTRACKVIEWANGLE", "MAXCROSSTRACKVIEWANGLE", "MEANCROSSTRACKVIEWANGLE",
                    "MINOFFNADIRVIEWANGLE", "MAXOFFNADIRVIEWANGLE", "MEANOFFNADIRVIEWANGLE",
                    "MINSUNAZ","MAXSUNAZ","MEANSUNAZ",
                    "MINSUNEL","MAXSUNEL","MEANSUNEL",
                    "MINSATAZ","MAXSATAZ","MEANSATAZ",
                    "MINSATEL","MAXSATEL","MEANSATEL",
                    "EARLIESTACQTIME","LATESTACQTIME",
                    "RADIOMETRICLEVEL", "RADIOMETRICENHANCEMENT"]

    for tag in metadata_tags:
        node = xmldoc.getElementsByTagName(tag)
        if len(node) > 0: metadata[tag] = node[0].firstChild.nodeValue

    if 'SATID' in metadata:
        if metadata['SATID'] == 'WV03':
            metadata['satellite'] = 'WorldView3'
            metadata['sensor'] = 'WorldView3'
            metadata["isotime"]=metadata["FIRSTLINETIME"]
            band_names=['COASTAL','BLUE','GREEN','YELLOW','RED','REDEDGE','NIR1','NIR2',
                        'SWIR1','SWIR2','SWIR3','SWIR4','SWIR5','SWIR6','SWIR7','SWIR8']
            band_indices=[1,2,3,4,5,6,7,8,
                          1,2,3,4,5,6,7,8]
            band_tag_names = ["BAND_C","BAND_B","BAND_G","BAND_Y","BAND_R","BAND_RE","BAND_N", "BAND_N2",
                              "BAND_S1", "BAND_S2", "BAND_S3", "BAND_S4", "BAND_S5", "BAND_S6", "BAND_S7", "BAND_S8"]
        if metadata['SATID'] == 'WV02':
            metadata['satellite'] = 'WorldView2'
            metadata['sensor'] = 'WorldView2'
            try:
                ## "L2A" data
                metadata["isotime"]=metadata["EARLIESTACQTIME"]
            except:
                ## "L1B" data
                metadata["isotime"]=metadata["FIRSTLINETIME"]
            band_names=['COASTAL','BLUE','GREEN','YELLOW','RED','REDEDGE','NIR1','NIR2']
            band_indices=[1,2,3,4,5,6,7,8]
            band_tag_names = ["BAND_C","BAND_B","BAND_G","BAND_Y","BAND_R","BAND_RE","BAND_N", "BAND_N2"]

        if metadata['SATID'] == 'QB02':
            metadata['satellite'] = 'QuickBird2'
            metadata['sensor'] = 'QuickBird2'
            try:
                ## "L2A" data
                metadata["isotime"]=metadata["EARLIESTACQTIME"]
            except:
                ## "L1B" data
                metadata["isotime"]=metadata["FIRSTLINETIME"]
            band_names=['Blue','Green','Red','NIR']
            band_indices=[1,2,3,4]
            band_tag_names = ["BAND_B","BAND_G","BAND_R","BAND_N"]

        if metadata['SATID'] == 'GE01':
            metadata['satellite'] = 'GeoEye1'
            metadata['sensor'] = 'GeoEye1'
            try:
                ## "L2A" data
                metadata["isotime"]=metadata["EARLIESTACQTIME"]
            except:
                ## "L1B" data
                metadata["isotime"]=metadata["FIRSTLINETIME"]
            band_names=['Blue','Green','Red','NIR']
            band_indices=[1,2,3,4]
            band_tag_names = ["BAND_B","BAND_G","BAND_R","BAND_N"]

    ## band tags to try and extract from metadata
    band_tags = ["ULLON","ULLAT","ULHAE",
                "URLON","URLAT","URHAE",
                "LRLON","LRLAT","LRHAE",
                "LLLON","LLLAT","LLHAE",
                "ABSCALFACTOR","EFFECTIVEBANDWIDTH","TDILEVEL"]

    ## read band information of spatial extent
    band_values={}
    band_index = 1
    for b,band_tag in enumerate(band_tag_names):
        #if band_tag == 'BAND_S1': band_index = 1 # reset counter for SWIR bands
        band_data = {'name':band_names[b], 'index':band_indices[b]}
        ## there are two tags in WV3 metadata
        for t in xmldoc.getElementsByTagName(band_tag):
            for tag in band_tags:
                    node = t.getElementsByTagName(tag)
                    if len(node) > 0:
                        band_data[tag]=float(node[0].firstChild.nodeValue)
        if len(band_data)>2: ## keep band only if in metadata
            ## change band index to the current index
            band_data['index'] = band_index
            band_index += 1
            ## copy over band data to dict
            band_values[band_tag]=band_data
    metadata['BAND_INFO'] = band_values

    ## get tile information
    tile_tags = ["FILENAME",
                "ULCOLOFFSET","ULROWOFFSET","URCOLOFFSET","URROWOFFSET",
                "LRCOLOFFSET","LRROWOFFSET","LLCOLOFFSET","LLROWOFFSET",
                "ULLON","ULLAT","URLON","URLAT",
                "LRLON","LRLAT","LLLON","LLLAT",
                "ULX","ULY","URX","URY",
                "LRX","LRY","LLX","LLY"]
    tile_values=[]
    for t in xmldoc.getElementsByTagName('TILE'):
        tile = {}
        for tag in tile_tags:
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    if tag == "FILENAME": val=node[0].firstChild.nodeValue
                    else: val=float(node[0].firstChild.nodeValue)
                    tile[tag]=val
        tile_values.append(tile)
    metadata['TILE_INFO'] = tile_values

    return(metadata)
