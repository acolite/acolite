## def copernicus_dem_download
## downloads Copernicus DEM TIFs or COGs from PRISM/AWS
## written by Quinten Vanhellemont, RBINS
## 2022-07-05
## modifications: 2023-09-12 (QV) added new Copernicus PRISM download options

def copernicus_dem_download(lon, lat, source = 'copernicus30'):
    import os, math, tarfile
    import acolite as ac
    local_dir = ac.config['copernicus_dem_dir']

    if source == 'copernicus30':
        resolution = 30
        base = 'Copernicus_DSM_COG_10'
        base_local = '{}'.format(base)
        url = 'https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com'
        server = 'amazon'
        ext = '_DEM'
    elif source == 'copernicus90':
        resolution = 90
        base = 'Copernicus_DSM_COG_30'
        base_local = '{}'.format(base)
        url = 'https://copernicus-dem-90m.s3.eu-central-1.amazonaws.com'
        server = 'amazon'
        ext = '_DEM'
    elif source == 'COP-DEM_GLO-30-DGED__2021_1':
        resolution = 30
        base = 'Copernicus_DSM_10'
        base_local = '{}'.format(source)
        url = 'https://prism-dem-open.copernicus.eu/pd-desk-open-access/prismDownload'
        server = 'copernicus'
        ext = ''
    elif source == 'COP-DEM_GLO-30-DGED__2022_1':
        resolution = 30
        base = 'Copernicus_DSM_10'
        base_local = '{}'.format(source)
        url = 'https://prism-dem-open.copernicus.eu/pd-desk-open-access/prismDownload'
        server = 'copernicus'
        ext = ''
    elif source == 'COP-DEM_GLO-90-DGED__2021_1':
        resolution = 90
        base = 'Copernicus_DSM_30'
        base_local = '{}'.format(source)
        url = 'https://prism-dem-open.copernicus.eu/pd-desk-open-access/prismDownload'
        server = 'copernicus'
        ext = ''
    elif source == 'COP-DEM_GLO-90-DGED__2022_1':
        resolution = 90
        base = 'Copernicus_DSM_30'
        base_local = '{}'.format(source)
        url = 'https://prism-dem-open.copernicus.eu/pd-desk-open-access/prismDownload'
        server = 'copernicus'
        ext = ''
    else:
        print('Copernicus DSM with identifier {} not configured.'.format(source))
        return()

    ## get tileList for this DSM
    tilelist = '{}/{}/tileList.txt'.format(local_dir, base_local)
    if (not os.path.exists(tilelist)):
        if (server == 'amazon'):
            download_url = '{}/tileList.txt'.format(url)
            ac.shared.download_file(download_url, tilelist)
        if (server == 'copernicus'):
            print('Getting Copernicus DEM file lists from PRISM. This can take a few minutes.')
            ac.dem.copernicus_dem_filelist()

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

    tile = '{}_{}{}_00_{}{}_00{}'.format(base, lat_char, str(lat_int).zfill(2), lon_char, str(lon_int).zfill(3), ext)

    ## test if tile exists
    if tile in tiles:
        if server == 'amazon':
            download_url = '{}/{}/{}.tif'.format(url, tile, tile)
            local_file = '{}/{}/{}'.format(local_dir, base_local, os.path.basename(download_url))
            if not os.path.exists(local_file):
                ac.shared.download_file(download_url, local_file)
        if server == 'copernicus':
            download_url = '{}/{}/{}.tar'.format(url, source, tile)
            local_file = '{}/{}/{}_DEM.tif'.format(local_dir, base_local, tile)
            if not os.path.exists(local_file):
                local_tar = '{}/{}/{}.tar'.format(local_dir, base_local, tile)
                if not os.path.exists(local_tar): ac.shared.download_file(download_url, local_tar)

                ## open tarfile to get DEM
                if os.path.exists(local_tar):
                    with tarfile.open(local_tar) as f:
                        files = f.getnames()
                        for fc in files:
                            if '{}/DEM/{}_DEM.tif'.format(tile,tile) in fc:
                                print('Extracting {}'.format(fc))
                                ## get member for this file
                                member = f.getmember(fc)
                                ## remove subdirectories
                                member.name = os.path.basename(member.name)
                                ## extract to target directory
                                f.extract(member, '{}/{}'.format(local_dir, base_local))
                ## delete tar file
                if os.path.exists(local_file) & os.path.exists(local_tar): os.remove(local_tar)
        return(local_file)
    else:
        print('Tile {} not in tilelist'.format(tile))
        return()
