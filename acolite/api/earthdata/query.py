## def query
## sets up queries for EarthData, downloads if asked
## written by Quinten Vanhellemont, RBINS
## 2024-04-28
## modifications:
##

def query(sensor, lon, lat, start_date = None, end_date = None, api = 'atom', verbosity = 5,
          download = False, local_directory = None, override = False,
          dataset = None, datacenter = None, collection_id = None,
          granule_id = None, ## search on granule name to be implemented "producergranuleid"
          filter_time = True, filter_time_range = [11, 14]): ## time filter for viirs to be implemented

    import os
    import acolite as ac

    if (download) & (local_directory is None):
        print('Please specify a local_directory for downloading scenes.')
        return

    apis = ['atom', 'json']
    if api not in apis:
        print('API {} not configured, use one of: {}'.format(api, ', '.join(apis)))
        return

    sensoru = sensor.upper()
    ## user can provide dataset and collection or use preconfigured sensor
    if (dataset is None) & (datacenter is None) & (collection_id is None):
        ## configure PACE
        if sensoru in ['PACE', 'OCI', 'PACE_OCI']:
            dataset = 'PACE_OCI_L1B_SCI'
            datacenter = 'OB_CLOUD'
            collection_id = 'C2804798373-OB_CLOUD' ## L1B
        elif sensoru in ['ECOSTRESS', 'ISS_ECOSTRESS']:
            collection_id = []
            collection_id.append('C1545228916-LPDAAC_ECS') ## ECO1BMAPRAD
            collection_id.append('C1534584923-LPDAAC_ECS') ## ECO1BGEO
            api = 'json'
        elif 'VIIRS' in sensoru:
            ## track which resolution files to get
            mod, img = True, True
            if 'MOD' in sensoru: img = False ## only MOD
            if 'IMG' in sensoru: mod = False ## only IMG

            ## if specific sensor is requested
            if 'SUOMI-NPP' in sensoru: sensor_ids = ['VNP']
            elif 'JPSS1' in sensoru: sensor_ids = ['VJ1']
            elif 'JPSS2' in sensoru: sensor_ids = ['VJ2']
            else: sensor_ids = ['VNP', 'VJ1', 'VJ2']

            dataset = []
            for sid in sensor_ids:
                if mod: dataset += ['{}02MOD'.format(sid), '{}03MOD'.format(sid)]
                if img: dataset += ['{}02IMG'.format(sid), '{}03IMG'.format(sid)]
            datacenter = ['LAADS'] * len(dataset)
            print(dataset)
            api = 'atom'
        else:
            print('Sensor {} not configured.'.format(sensor))
            return

    ## convert to list since we need to iterate for sensors requiring multiple files
    if type(dataset) is not list: dataset = [dataset]
    if type(datacenter) is not list: datacenter = [datacenter]
    if type(collection_id) is not list: collection_id = [collection_id]
    n = max((len(dataset), len(collection_id)))

    if verbosity > 3:
        print('Using dataset {}, datacenter {}, collection_id {}'.format(dataset, datacenter, collection_id))

    ## find urls and scenes for each requested dataset/collection_id
    urls, files = [], []
    for di in range(n):
        ## construct query url
        if api == 'atom':
            query_url = ac.api.earthdata.url_atom(dataset[di], lat, lon, start_date = start_date,
                                                    end_date = end_date, datacenter = datacenter[di])
        elif api == 'json':
            query_url = ac.api.earthdata.url_json(collection_id[di], lat, lon,
                                                    start_date = start_date, end_date = end_date)
        if verbosity > 3: print(query_url)

        ## get scenes and urls from query url
        urls_, files_ = ac.api.earthdata.query_extract_scenes(query_url, verbosity = verbosity, link_base = None)
        if verbosity > 3: print(urls_, files_)
        urls+=urls_
        files+=files_
    ## end run through different queries

    ## return urls and files when not downloading
    if not download: return(urls, files)

    ## download files if we do not have a local copy
    local_files = []
    for i, url in enumerate(urls):
        scene = os.path.basename(url)
        if scene != files[i]:
            if verbosity > 1:
                print('URL does not match scene')
                print(url)
                print(files[i], scene)
            continue

        ## local file path
        file = '{}/{}'.format(local_directory, os.path.basename(url))
        if (not os.path.exists(file)) | (override):
            if verbosity > 1: print('Downloading {}'.format(url))
            ac.shared.download_file(url, file)
        else:
            if verbosity > 1: print('We have {}'.format(url))

        if verbosity > 1: print('Local file {}'.format(file))
        local_files.append(file)
    return(local_files)
