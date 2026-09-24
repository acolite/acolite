## def query
## finds Planet data for a given date range and region of interest
## written by Quinten Vanhellemont, RBINS
## 2020-10-13
## modifications: 2022-01-11 (QV) update for SuperDove data
##                2026-09-22 (QV) integrated in acolite.api.planet

def query(geojson_geometry, date_start, date_end = None,
                                output = None,
                                cloud_filter = True, max_cloud = 0.5, min_cover = 0.5,
                                item_type = "PSScene", #"PSScene5Band", #"PSOrthoTile", #"PSScene5Band" PSScene4Band", #"REOrthoTile"
                                search_url = 'https://api.planet.com/data/v1/quick-search',
                                orders_url = 'https://api.planet.com/compute/ops/orders/v2',
                                overwrite = False, filter_ids = None):

    import json, os, pathlib, requests
    from requests.auth import HTTPBasicAuth
    from shapely.geometry import Polygon
    import acolite as ac

    # API Key stored as an env variable
    if 'PL_API_KEY' in os.environ:
        PL_API_KEY = os.getenv('PL_API_KEY')
    else:
        print('PL_API_KEY not in environment. Please add your Planet API key as PL_API_KEY')
        return

    # set up requests to work with api
    auth = HTTPBasicAuth(PL_API_KEY, '')
    headers = {'content-type': 'application/json'}

    ## set up end date
    if date_end is None: date_end = date_start ## one millisecond before midnight

    # get images that overlap with our AOI
    geometry_filter = {
      "type": "GeometryFilter",
      "field_name": "geometry",
      "config": geojson_geometry
    }

    ## roi Polygon
    roi = Polygon(geojson_geometry['coordinates'][0])

    # get images acquired within a date range
    date_range_filter = {
      "type": "DateRangeFilter",
      "field_name": "acquired",
      "config": {
        "gte": "{}T00:00:00.000Z".format(date_start),
        "lte": "{}T23:59:59.999Z".format(date_end)
      }
    }

    # combine filters
    combined_filter = {
      "type": "AndFilter",
      "config": [geometry_filter, date_range_filter]}

    ## cloud filter
    if cloud_filter:
        combined_filter['config'].append({"type": "RangeFilter", "field_name": "cloud_cover",
                                          "config": {"lte": max_cloud}})

    # API request object
    search_request = {
      "item_types": [item_type],
      "filter": combined_filter
    }

    #if True:
    #    search_request['instrument_type'] = 'PSD.SD'

    # send POST request
    search_result = requests.post(search_url,auth=HTTPBasicAuth(PL_API_KEY, ''), json=search_request)
    files = []
    features_ = search_result.json()['features']

    ## filter on id if provided
    features = []
    for f in features_:
        if type(filter_ids) is list:
            if f['id'] in filter_ids: continue
        features.append(f)
    print('Found {} features'.format(len(features)), end='\n')

    ## return retrieved features
    if min_cover is not None:
        covered_features = []
        for f in features:
            scene_extent = Polygon(f['geometry']['coordinates'][0])
            x = roi.intersection(scene_extent)
            fr = x.area/roi.area
            if fr>= min_cover:
                covered_features.append(f)

        print('Found {} features with cloud cover < {:.2f}, covering ROI > {:.2f}'.format(len(covered_features), max_cloud, min_cover), end='\n')
        return(covered_features)
    else:
        return(features)
