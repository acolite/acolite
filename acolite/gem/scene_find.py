## def scene_find
## find (and download) scenes for given position and date range
## written by Quinten Vanhellemont, RBINS
## 2021-03-11
## modifications:

def scene_find(st_lon, st_lat, sdate, edate=None,
                sources = ['Landsat 5', 'Landsat 7','Landsat 8','Sentinel-2'],
                output_base=None, download=True, verbosity=0):

    import acolite as ac

    iml = ac.gem.extract(st_lon, st_lat, sdate, edate=edate, return_iml=True, sources=sources,
                         override=False,st_name = None, output=None, verbosity=verbosity)

    dr = '{}'.format(sdate) if edate is None else '{}-{}'.format(sdate, edate)
    if verbosity > 0: print('Found {} scene{} for {}N {}E ({})'.format(len(iml), '' if len(iml) == 1 else 's', st_lat, st_lon, dr))

    local_scenes = []
    for im in iml:
        ret = ac.gem.scene_download(im, output_base=output_base, download=download, verbosity=verbosity)
        local_scenes.append(ret)

    return(local_scenes)
