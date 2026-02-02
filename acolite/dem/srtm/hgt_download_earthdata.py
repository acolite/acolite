## def srtm.hgt_download_earthdata
## downloads SRTM granule from EarthData
## written by Quinten Vanhellemont, RBINS
## 2026-02-02
## modifications:

def hgt_download_earthdata(granule, source = None, verbosity = 6, download = True, override = False, output = None):
    import requests, time, json, os
    import acolite as ac

    #datacenter = 'LPCLOUD'
    #api = 'json'
    if granule.endswith('.zip'):
        granule = granule[0:granule.find('.zip')]
    source = granule.split('.')[1]

    ## find collection ids
    if source == 'SRTMGL3':
        collection_id = 'C2763266377-LPCLOUD'
        dataset = 'NASA Shuttle Radar Topography Mission Global 3 arc second V003'
    elif source == 'SRTMGL3S':
        collection_id = 'C2763268442-LPCLOUD'
        dataset = 'NASA Shuttle Radar Topography Mission Global 3 arc second sub-sampled V003'
    elif source == 'SRTMGL1':
        collection_id = 'C2763266360-LPCLOUD'
        dataset = 'NASA Shuttle Radar Topography Mission Global 1 arc second V003'
    else:
        print('Source={} not configured'.format(source))
        return

    if output is None:
        hgt_dir = '{}/{}'.format(ac.config['hgt_dir'], source)
    else:
        hgt_dir = '{}'.format(output)

    ## set up query
    query_base = 'https://cmr.earthdata.nasa.gov:443/search/granules.json'
    query_url = '{}?'.format(query_base)
    query_url += '&echo_collection_id={}'.format(collection_id)
    query_url += '&granuleUr={}'.format(granule)

    ## post query
    if verbosity > 0: print('Posting query {}'.format(query_url))
    response = requests.get(query_url)
    if response.ok:
        dct = json.loads(response.text)
    else:
        print('Could not load EarthData query response.')
        return

    ## run trough query entries
    urls = []
    for entry in dct['feed']['entry']:
        for link in entry['links']:
            if ('href' in link.keys()) & ('title' in link.keys()):
                if granule + '.zip' in link['title']:
                    urls.append(link['href'])

    ## return urls and files when not downloading
    if not download: return(urls)

    ## download files if we do not have a local copy
    local_files = []
    for i, url in enumerate(urls):
        ## local file path
        file = '{}/{}'.format(hgt_dir, os.path.basename(url))
        if (not os.path.exists(file)) | (override):
            if verbosity > 1: print('Downloading {}'.format(url))
            ac.shared.download_file(url, file)
        else:
            if verbosity > 1: print('We have {}'.format(url))

        if verbosity > 1: print('Local file {}'.format(file))
        local_files.append(file)
    return(local_files)
