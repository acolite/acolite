## def polylakes
## downloads and extracts a global lake polygon shapefile, returns local path
##
## World Lakes is provided by Daniel Siegel at the ArcGIS Living Atlas Team
## under CCv40 license https://creativecommons.org/licenses/by/4.0/
## link: https://hub.arcgis.com/datasets/0abb136c398942e080f736c8eb09f5c4_0/about
##
## HydroLAKES is provided by HydroSHEDS
## link: https://www.hydrosheds.org/page/hydrolakes
## reference: Messager et al. 2016 doi: 10.1038/ncomms13603
##
## function written by Quinten Vanhellemont, RBINS
## 2022-01-01
## modifications: 2022-01-01 (QV) renamed from worldlakes, added hydrolakes

def polylakes(database = 'worldlakes', remove_zip = False):
    import acolite as ac
    import os, zipfile

    if database.lower() not in ['worldlakes', 'hydrolakes']:
        print('Polylake database {} not recognised, using worldlakes.')
        database = 'worldlakes'

    ## local worldlakes files
    if database.lower() == 'worldlakes':
        url = 'https://opendata.arcgis.com/api/v3/datasets/0abb136c398942e080f736c8eb09f5c4_0/downloads/data?format=shp&spatialRefId=4326'
        local_zip = ac.config['path']+'/external/{}'.format('World_Lakes-shp.zip')
        local_file = ac.config['path']+'/external/{}'.format('World_Lakes-shp/World_Lakes.shp')
        local_dir = ac.config['path']+'/external/{}'.format('World_Lakes-shp')
    elif database.lower() == 'hydrolakes':
        url = 'https://97dc600d3ccc765f840c-d5a4231de41cd7a15e06ac00b0bcc552.ssl.cf5.rackcdn.com/HydroLAKES_polys_v10_shp.zip'
        local_zip = ac.config['path']+'/external/{}'.format('HydroLAKES_polys_v10_shp.zip')
        local_file = ac.config['path']+'/external/{}'.format('HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp')
        local_dir = ac.config['path']+'/external/{}'.format('HydroLAKES_polys_v10_shp')

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
                z.extractall(local_dir)
            if remove_zip: os.remove(local_zip)

    ## return local path
    if os.path.exists(local_file):
        return(local_file)
    else:
        return(None)
