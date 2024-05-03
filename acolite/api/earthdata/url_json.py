## def url_json
## get query url using the json search
## written by Quinten Vanhellemont, RBINS
## 2023-11-22
## modifications: 2024-04-15 (QV) renamed from query
##                2024-04-28 (QV) added as acolite function, renamed from query_json, added scene search option

def url_json(collection_id, lat, lon, scene = None, start_date = None, end_date = None,
          query_base = 'https://cmr.earthdata.nasa.gov:443/search/granules.json',
          ): # nresults = 1000, radius = 1000 not used at the moment

    import datetime, dateutil.parser

    if scene is not None: ## search for scene
        query_url = '{}?'.format(query_base)
        query_url += '&echo_collection_id={}'.format(collection_id)
        query_url += '&producerGranuleId={}'.format(scene)
    else:  ## search for location and datetime
        query_url = '{}?'.format(query_base)
        query_url += '&echo_collection_id={}'.format(collection_id)
        query_url += '&point={},{}'.format(lon,lat)

        temporal = ['','','','']
        if start_date is not None:
            sdate = dateutil.parser.parse(start_date)
            temporal[0] = '{}Z'.format(sdate.isoformat())
            #query_url += '&startTime={}Z'.format(sdate.isoformat())

        if end_date is not None:
            edate = dateutil.parser.parse(end_date)
            add_day = 0
            if start_date is not None:
                if edate == sdate: add_day = 1 ## add one day if start = end date
            if 'T' not in end_date: add_day = 1 ## add one day if no time is given
            edate+=datetime.timedelta(days=add_day)
            temporal[1] = '{}Z'.format(edate.isoformat())
        query_url += '&temporal={}'.format(','.join([str(v) for v in temporal]))

    query_url=query_url.replace(' ', '%20')
    return(query_url)
