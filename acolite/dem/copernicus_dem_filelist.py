## def copernicus_prism_dem_filelist
## gets DEM file lists from PRISM server
## written by Quinten Vanhellemont, RBINS
## 2023-09-12
## modifications:

def copernicus_dem_filelist(publicDemURLs = 'https://prism-dem-open.copernicus.eu/pd-desk-open-access/publicDemURLs'):

    import os, requests
    import acolite as ac
    local_dir = ac.config['copernicus_dem_dir']


    ret = requests.get(publicDemURLs, headers={"accept":"json"})
    datasets = [r['datasetId'] for r in ret.json()]

    for dataset in datasets:
        base = dataset.replace('/','__')
        tilelist = '{}/{}/tileList.txt'.format(local_dir, base)
        if not os.path.exists(os.path.dirname(tilelist)):
            os.makedirs(os.path.dirname(tilelist))

        if not os.path.exists(tilelist):
            ret = requests.get(publicDemURLs + '/' + base, headers={"accept":"json"})
            urls = [r['nativeDemUrl'] for r in ret.json()]

            ## find tiles
            tiles = []
            for url in urls:
                tile, ext = os.path.splitext(os.path.basename(url))
                tiles.append(tile)
            tiles.sort()

            ## write tilelist
            with open(tilelist, 'w') as f:
                for tile in tiles:
                    f.write(tile)
                    f.write('\n')
            print('Wrote {}'.format(tilelist))
