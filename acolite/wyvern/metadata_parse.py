## def metadata_parse
## gets required info from Wyvern json file
## written by Quinten Vanhellemont, RBINS
## 2025-03-04
## modifications:

def metadata_parse(metafile):
    import json
    import acolite as ac
    import numpy as np

    ## read meta file
    with open(metafile, 'r', encoding='utf-8') as f:
        meta = json.load(f)

    gatts = {}
    gatts['sensor'] = '{}_{}'.format(meta['properties']['processing:facility'], meta['properties']['platform'])
    gatts['isodate'] = meta['properties']['datetime']
    gatts['resolution'] = meta['properties']['gsd']

    gatts['sza'] = 90 - meta['properties']['view:sun_elevation']
    gatts['saa'] = meta['properties']['view:sun_azimuth']
    gatts['vza'] = meta['properties']['view:incidence_angle']
    gatts['vaa'] = meta['properties']['view:azimuth']
    gatts['raa'] = abs(gatts['saa'] - gatts['vaa'])
    while gatts['raa'] > 180: gatts['raa'] = np.abs(gatts['raa']-360)

    return(meta, gatts)
