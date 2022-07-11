## def copernicus_dem_download
## downloads AWS Copernicus COGs
## written by Quinten Vanhellemont, RBINS
## 2022-07-05

def copernicus_dem_download(lon, lat, source = 'copernicus30'):
    import os, math
    import acolite as ac
    local_dir = ac.config['copernicus_dem_dir']
    
    if source == 'copernicus30':
        resolution = 30
        base = 'Copernicus_DSM_COG_10'
        url = 'https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com'
    elif source == 'copernicus90':
        resolution = 90
        base = 'Copernicus_DSM_COG_30'
        url = 'https://copernicus-dem-90m.s3.eu-central-1.amazonaws.com'
    else:
        print('Copernicus DSM with resolution {} not configured.'.format(resolution))
        return()

    ## get tileList for this DSM
    tilelist = '{}/{}/tileList.txt'.format(local_dir, base)
    if not os.path.exists(tilelist):
        download_url = '{}/tileList.txt'.format(url)
        ac.shared.download_file(download_url, tilelist)
    ## read tiles
    tiles = ac.acolite.settings.read_list(tilelist)

    ## find tile needed for current position
    lat_int = abs(int(math.floor(lat)))
    lon_int = abs(int(math.floor(lon)))
    if lat < 0:
        lat_char = 'S'
    else:
        lat_char = 'N'
    if lon < 0:
        lon_char = 'W'
    else:
        lon_char = 'E'
    tile = '{}_{}{}_00_{}{}_00_DEM'.format(base, lat_char, str(lat_int).zfill(2), lon_char, str(lon_int).zfill(3))

    ## test if tile exists
    if tile in tiles:
        download_url = '{}/{}/{}.tif'.format(url, tile, tile)
        local_file = '{}/{}/{}'.format(local_dir, base, os.path.basename(download_url))
        if not os.path.exists(local_file):
            ac.shared.download_file(download_url, local_file)
        return(local_file)
    else:
        print('Tile {} not in tilelist'.format(tile))
        return()
