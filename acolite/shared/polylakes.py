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
##                2024-02-29 (QV) added external dir config
##                                added gshhg (L1 for land/water mask)

def polylakes(database = 'worldlakes', remove_zip = False):
    import acolite as ac
    import os, zipfile

    if database.lower() not in ['worldlakes', 'hydrolakes', 'gshhg']:
        print('Polylake database {} not recognised, using worldlakes.'.format(database))
        database = 'worldlakes'

    ## local worldlakes files
    if database.lower() == 'worldlakes':
        url = 'https://opendata.arcgis.com/api/v3/datasets/0abb136c398942e080f736c8eb09f5c4_0/downloads/data?format=shp&spatialRefId=4326'
        local_zip = '{}/{}'.format(ac.config['directory']['external'], 'World_Lakes-shp.zip')
        local_file = '{}/{}'.format(ac.config['directory']['external'], 'World_Lakes-shp/World_Lakes.shp')
        local_dir = '{}/{}'.format(ac.config['directory']['external'], 'World_Lakes-shp')
    elif database.lower() == 'hydrolakes':
        url = 'https://97dc600d3ccc765f840c-d5a4231de41cd7a15e06ac00b0bcc552.ssl.cf5.rackcdn.com/HydroLAKES_polys_v10_shp.zip'
        local_zip = '{}/{}'.format(ac.config['directory']['external'], 'HydroLAKES_polys_v10_shp.zip')
        local_file = '{}/{}'.format(ac.config['directory']['external'], 'HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp')
        local_dir = '{}/{}'.format(ac.config['directory']['external'], 'HydroLAKES_polys_v10_shp')
    elif database.lower() == 'gshhg':
        url = 'https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip'
        local_zip = '{}/{}'.format(ac.config['directory']['external'], 'gshhg-shp-2.3.7.zip')
        local_file = '{}/{}'.format(ac.config['directory']['external'], 'gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp')
        local_dir = '{}/{}'.format(ac.config['directory']['external'], 'gshhg-shp-2.3.7')

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
