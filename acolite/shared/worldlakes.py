## def worldlakes
## downloads and extracts the world lakes shapefile, returns local path
##
## World Lakes is provided by Daniel Siegel at the ArcGIS Living Atlas Team
## under CCv40 license https://creativecommons.org/licenses/by/4.0/
## link: https://hub.arcgis.com/datasets/0abb136c398942e080f736c8eb09f5c4_0/about
##
## function written by Quinten Vanhellemont, RBINS
## 2022-01-01
## modifications:

def worldlakes(url = 'https://opendata.arcgis.com/api/v3/datasets/0abb136c398942e080f736c8eb09f5c4_0/downloads/data?format=shp&spatialRefId=4326',
               remove_zip = False):
    import acolite as ac
    import os, zipfile

    ## local worldlakes files
    local_zip = ac.config['path']+'/external/{}'.format('World_Lakes-shp.zip')
    local_file = ac.config['path']+'/external/{}'.format('World_Lakes-shp/World_Lakes.shp')

    ## download/extract
    if not os.path.exists(local_file):
        ## download
        if not os.path.exists(local_zip):
            print('Downloading {}'.format(local_zip))
            ret = ac.shared.download_file(url, local_zip)

        ## extract
        if os.path.exists(local_zip):
            print('Extracting {}'.format(local_zip))
            with zipfile.ZipFile(local_zip, 'r') as z:
                z.extractall('{}/external/{}/'.format(ac.config['path'], 'World_Lakes-shp'))
            if remove_zip: os.remove(local_zip)

    ## return local path
    if os.path.exists(local_file):
        return(local_file)
    else:
        return(None)
