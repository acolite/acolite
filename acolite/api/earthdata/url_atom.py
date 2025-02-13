## def url_atom
## get query url using the atom search
## written by Quinten Vanhellemont, RBINS
## 2023-10-26
## modifications: 2024-04-15 (QV) renamed from query, added datacenter as keyword
##                2024-04-28 (QV) added as acolite function, renamed from query_atom
##                2025-02-13 (QV) added bounding box option

def url_atom(dataset, lat, lon, limit = None, scene = None, start_date = None, end_date = None, datacenter = 'LAADS',
          query_base='https://cmr.earthdata.nasa.gov/opensearch/granules.atom', nresults = 1000, radius = 1000):

    import datetime, dateutil.parser
    import numpy as np

    query_url = '{}?'.format(query_base)
    query_url += '&shortName={}&dataCenter={}'.format(dataset, datacenter)

    if scene is not None: ## not yet working as intended
        query_url += '&producerGranuleId={}'.format(scene)
    else:
        ## one point
        if (np.atleast_1d(lat).shape[0] == 1) & (np.atleast_1d(lon).shape[0] == 1):
            query_url += '&lat={}'.format(lat)
            query_url += '&lon={}'.format(lon)
            query_url += '&radius={}'.format(radius)
        ## bounding box
        elif (np.atleast_1d(lat).shape[0] == 2) & (np.atleast_1d(lon).shape[0] == 2):
            query_url += '&boundingBox={:.5f},{:.5f},{:.5f},{:.5f}'.format(lon[0],lat[0],lon[1],lat[1])
        else:
            if (limit is not None):
                if np.atleast_1d(limit).shape[0] == 4:
                    query_url += '&boundingBox={:.5f},{:.5f},{:.5f},{:.5f}'.format(limit[1],limit[0],limit[3],limit[2])
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
