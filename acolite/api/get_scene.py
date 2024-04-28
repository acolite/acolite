## def get_scene
## try to get scene from a download API
## written by Quinten Vanhellemont, RBINS
## 2024-04-28
## modifications: 2024-04-28 (QV) split off from acolite.inputfile_test

def get_scene(scene, download_directory = None):
    import os
    import acolite as ac

    ## identify scene
    bn = os.path.basename(scene)
    if bn[0:3] in ['S2A', 'S2B', 'S3A', 'S3B']:
        download_source = 'CDSE'
    elif bn[0:4] in ['LC08', 'LO08', 'LT08', 'LC09', 'LO09', 'LT09', 'LT04', 'LT05', 'LE07']:
        download_source = 'EarthExplorer'
    elif 'ECOSTRESS' in bn:
        download_source = 'EarthExplorer'
    elif ('PACE_OCI' in bn):
        download_source = 'EarthData'
        sensor = 'PACE_OCI'
    elif (bn[0:3] in ['VNP', 'VJ1', 'VJ2']):
        download_source = 'EarthData'
        sensor = bn[0:3]
        print('VIIRS scene download not yet implemented')
        return
    else:
        print('Could not identify download source for scene {}'.format(file))
        return

    ##
    if ac.config['verbosity'] > 0: print('Attempting download of scene {} from {}.'.format(scene, download_source))
    if ac.config['verbosity'] > 0: print('Querying {}'.format(download_source))

    ## Copernicus Data Space Ecosystem
    if download_source == 'CDSE':
        urls, scenes = ac.api.cdse.query(scene=bn)
        if ac.config['verbosity'] > 0: print('Downloading from {}'.format(download_source))
        local_scenes = ac.api.cdse.download(urls, output = download_directory, verbosity = ac.config['verbosity'])

    ## EarthExplorer
    if download_source == 'EarthExplorer':
        entity_list, identifier_list, dataset_list = ac.api.earthexplorer.query(scene=bn)
        if ac.config['verbosity'] > 0: print('Downloading from {}'.format(download_source))
        local_scenes = ac.api.earthexplorer.download(entity_list, dataset_list, identifier_list,
                                                     output = download_directory, verbosity = ac.config['verbosity'])

    ## EarthData
    if download_source == 'EarthData':
        local_scenes = ac.api.earthdata.query(sensor, scene = bn, download = True,
                                              local_directory = download_directory, verbosity = ac.config['verbosity'])

    ## return local paths
    return(local_scenes)
