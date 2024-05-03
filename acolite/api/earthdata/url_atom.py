## def url_atom
## get query url using the atom search
## written by Quinten Vanhellemont, RBINS
## 2023-10-26
## modifications: 2024-04-15 (QV) renamed from query, added datacenter as keyword
##                2024-04-28 (QV) added as acolite function, renamed from query_atom

def url_atom(dataset, lat, lon, scene = None, start_date = None, end_date = None, datacenter = 'LAADS',
          query_base='https://cmr.earthdata.nasa.gov/opensearch/granules.atom', nresults = 1000, radius = 1000):

    import datetime, dateutil.parser

    query_url = '{}?'.format(query_base)
    query_url += '&shortName={}&dataCenter={}'.format(dataset, datacenter)

    if scene is not None: ## not yet working as intended
        query_url += '&producerGranuleId={}'.format(scene)
    else:
        query_url += '&lat={}'.format(lat)
        query_url += '&lon={}'.format(lon)
        query_url += '&radius={}'.format(radius)
        query_url += '&numberOfResults={}'.format(nresults)

        if start_date is not None:
            sdate = dateutil.parser.parse(start_date)
            query_url += '&startTime={}Z'.format(sdate.isoformat())

        if end_date is not None:
            edate = dateutil.parser.parse(end_date)
            add_day = 0
            if start_date is not None:
                if edate == sdate: add_day = 1 ## add one day if start = end date
            if 'T' not in end_date: add_day = 1 ## add one day if no time is given
            edate+=datetime.timedelta(days=add_day)
            query_url += '&endTime={}Z'.format(edate.isoformat())

    return(query_url)
