## def query
## create stac query for S2 zarr data
## written by Quinten Vanhellemont, RBINS
## 2025-11-17
## modifications:

def query(limit = None, start_date = None, end_date = None, nlimit = 100, scene = None, sensor = None,
            select_features = False, base_url = 'https://stac.core.eopf.eodc.eu/', collection = 'sentinel-2-l1c'):

    import json, requests, dateutil.parser, datetime
    import acolite as ac
    import numpy as np

    ## set up base query url
    query_url = '{}/collections/{}/items'.format(base_url, collection)

    ## add bounding box to query
    if limit is not None:
        query_url+='?bbox={},{},{},{}'.format(limit[1],limit[0],limit[3], limit[2])

    ## add datetime to query:
    date_query = ''
    if start_date is None:
        date_query += '/'
        sdate = None
    else:
        sdate = dateutil.parser.parse(start_date)
        if sdate.tzinfo is None: sdate = sdate.replace(tzinfo=datetime.timezone.utc)
        sdate = sdate.astimezone(datetime.timezone.utc)
        date_query += '{}/'.format(sdate.isoformat()[0:19]+'Z')

    if end_date is not None:
        edate = dateutil.parser.parse(end_date)
        if edate.tzinfo is None: edate = edate.replace(tzinfo=datetime.timezone.utc)
        edate = edate.astimezone(datetime.timezone.utc)
        if edate == sdate:
            edate += datetime.timedelta(days = 1)
        date_query += '{}'.format(edate.isoformat()[0:19]+'Z')

    if len(date_query) > 1:
        if '?' not in query_url:
            query_url = '?'+query_url
        else:
            query_url+='&'
        query_url+='datetime={}'.format(date_query)


    ## add max number of items
    query_url += '&limit={}'.format(nlimit)

    ## directly create url if scene is given
    if scene is not None:
        query_url = '{}/collections/{}/items/{}'.format(base_url, collection, scene)

    ## query features
    features = ac.api.stac.query_get(query_url)

    if select_features:
        selected_features = []
        for f in features:
            if sensor is not None:
                if not f['id'].startswith(sensor): continue
            zipped_product = None

            product = f['assets']['product']['href']
            if 'zipped_product' in f['assets']:
                zipped_product = f['assets']['zipped_product']['href']
            selected_features.append(f)

        sort = np.argsort([f['id'] for f in selected_features])
        features_ids = [selected_features[i]['id'] for i in sort]
        feature_urls = [selected_features[i]['assets']['product']['href'] for i in sort]
        return(features_ids, feature_urls)

    return(features)
