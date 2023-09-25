## def query
## queries EarthExplorer for landsat scenes
## written by Quinten Vanhellemont, RBINS
## 2023-09-19
## modifications: 2023-09-20 (QV) moved roi_wkt to function, added cloudCoverFilter

def query(scene = None, collection = 2, level = 1, landsat_type = None,
               start_date = None, end_date = None,  roi = None,
               verbosity = 1, cloud_cover = None,
               max_results = 1000, api_url = None, attributes = False):

    import os, requests, json, netrc
    import datetime, dateutil.parser
    from osgeo import ogr,osr,gdal
    import acolite as ac

    ## get api URL from config
    if api_url is None: api_url = ac.config['EARTHEXPLORER_api']

    ## get access token for query
    access_token = ac.earthexplorer.auth(api_url=api_url)
    if access_token is None:
        print('Could not get EarthExplorer access token!')
        return

    ## get landsat_type from scene
    if scene is not None:
        landsat_type = None
        if scene[3] in ['8','9']:
            landsat_type = 'ot'
        if scene[3] in ['4','5']:
            landsat_type = 'tm'
        if scene[3] in ['7']:
            landsat_type = 'etm'
        if landsat_type is None:
            print('Could not determine landsat_type based on scene name "{}"'.format(scene))
            return
        print('Assuming landsat_type "{}" based on scene name'.format(scene))

    ## set up dataset
    landsat_type_default = 'ot'
    if landsat_type is None:
        print('Assuming landsat_type "ot" (Landsat 8, Landsat 9 OLI/TIRS)')
        print('Options for landsat_type: "ot" (L8/L9),  "tm" (L4/L5), "etm" (L7)')
        landsat_type = '{}'.format(landsat_type_default)

    if landsat_type not in ['ot', 'tm', 'etm']:
        print('landsat_type "{}" not supported'.format(landsat_type))
        print('Using landsat_type "ot" (Landsat 8, Landsat 9 OLI/TIRS)')
        print('Options for landsat_type: "ot" (L8/L9),  "tm" (L4/L5), "etm" (L7)')
        landsat_type = '{}'.format(landsat_type_default)

    if level not in [1,2]:
        print('level {} not supported, using level 1'.format(level))
        level = 1

    if collection not in [2]:
        print('collection {} not supported, using collection 2'.format(collection))
        collection = 2
    ## compose dataset
    dataset = 'landsat_{}_c{}_l{}'.format(landsat_type, collection, level)

    print('Using dataset {}'.format(dataset))

    ## determine WKT from provided ROI
    wkt = None
    if (roi is not None): wkt = ac.shared.roi_wkt(roi)
    if wkt is not None:
        if verbosity > 1: print('Using WKT for query: {}'.format(wkt))

    ## empty lists for results
    ## this could be neater, but works
    entity_list, dataset_list = [], []
    metadata_list = []

    ## set up session
    session = requests.Session()
    session.headers["X-Auth-Token"] = f"{access_token}"

    ## query for scene
    if scene is not None:
        ## find old style Landsat scene IDs
        data = {"listId": 'acolite_query', "datasetName": dataset, "idField": "displayId", "entityId": scene}
        response = session.post(api_url+'/scene-list-add', data = json.dumps(data))

        ## get scene list
        response = session.post(api_url+'/scene-list-get', data = json.dumps({"listId": 'acolite_query'}))

        ## parse scene list, get metadata
        for s in response.json()['data']:
            if verbosity > 2: print(s['entityId'])
            entity_list.append(s['entityId'])
            dataset_list.append(s['datasetName'])

        ## delete scene list
        response = session.post(api_url+'/scene-list-remove', data = json.dumps({"listId": 'acolite_query'}))

        ## find scene metadata - optional
        for si, s in enumerate(entity_list):
            data = {"datasetName": dataset_list[si], "entityId": entity_list[si], "metadataType": "full"}
            response = session.post(api_url+'/scene-metadata', data = json.dumps(data))
            metadata_list.append(response.json()['data']['metadata'])

    else:
        query = {}
        #if (lat is not None) and (lon is not None):
        #    coord = {"longitude": lon, "latitude": lat}
        #    query['spatialFilter'] = {'filterType': 'mbr', 'lowerLeft': coord, 'upperRight': coord}

        if (wkt is not None):
            if verbosity > 1: print('Finding boundary of WKT: {}'.format(wkt))

            geom = ogr.CreateGeometryFromWkt(wkt)
            if 'POINT' in wkt:
                points = geom.GetPoints()
                coordLL = {"longitude": points[0][0], "latitude": points[0][1]}
                coordUR = {"longitude": points[0][0], "latitude": points[0][1]}
            else:
                points = geom.GetBoundary().GetPoints()
                coordLL = {"longitude": points[0][0], "latitude": points[0][1]}
                coordUR = {"longitude": points[2][0], "latitude": points[2][1]}
            if verbosity > 1: print('Using boundary LL: {}N {}E'.format(coordLL['latitude'], coordLL['longitude']))
            if verbosity > 1: print('Using boundary UR: {}N {}E'.format(coordUR['latitude'], coordUR['longitude']))

            geom = None
            query['spatialFilter'] = {'filterType': 'mbr', 'lowerLeft': coordLL, 'upperRight': coordUR}

        if (start_date is not None) | (end_date is not None):
            query['acquisitionFilter'] = {}
            if start_date is not None:
                sdate = dateutil.parser.parse(start_date)
                query['acquisitionFilter']['start'] = sdate.isoformat()

            if end_date is not None:
                edate = dateutil.parser.parse(end_date)
                edate += datetime.timedelta(days=1) ## add one day to include end date data
                query['acquisitionFilter']['end'] = edate.isoformat()

        if (cloud_cover is not None):
            query['cloudCoverFilter'] = {'max': cloud_cover, 'min': 0, 'includeUnknown': False}

        data = {"datasetName": dataset, "sceneFilter": query, "maxResults": max_results, "metadataType": "full"}
        if verbosity > 3: print(data)
        response = session.post(api_url+'/scene-search', data = json.dumps(data))

        ## parse results metadata
        if (response.ok):
            for s in response.json()['data']['results']:
                if verbosity > 2: print(s['entityId'])
                entity_list.append(s['entityId'])
                dataset_list.append(dataset)
                metadata_list.append(s['metadata'])
        else:
            print("Error in posting query.")

    ## get product identifier from metadata
    identifier_list = []
    for sc in metadata_list:
        for m in sc:
             if 'Product Identifier' in m['fieldName']:
                identifier_list.append(m['value'])
    if verbosity > 0: print('Found {} total scenes'.format(len(entity_list)))
    if attributes: return(entity_list, identifier_list, dataset_list, metadata_list)
    return(entity_list, identifier_list, dataset_list)