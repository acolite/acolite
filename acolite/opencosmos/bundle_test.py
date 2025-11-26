## def bundle_test
## idendifies sensor and finds json and tiff files for OpenCosmos data
## currently Menut and Platero
## written by Quinten Vanhellemont, RBINS
## 2025-03-14
## modifications: 2025-04-28 (QV) update indexing to get band name
##                2025-11-26 (QV) added HAMMER

def bundle_test(bundle):
    import acolite as ac
    import os, glob

    if os.path.isdir(bundle):
        files = glob.glob('{}/*'.format(bundle))
        files.sort()
    elif os.path.isfile(bundle):
        files = glob.glob('{}/*'.format(os.path.dirname(bundle)))
        files.sort()
    else:
        return

    basenames = [os.path.basename(f) for f in files]

    ## find tiff files and identify available bands
    images = [b for b in basenames if b.endswith('.tiff')]
    bands = [i.split('_')[-2] for i in images]
    sensor = None

    ## identify sensor - should be improved
    menut_bands = ['B', 'G', 'R', 'RE1', 'RE2', 'RE3', 'NIR']
    platero_bands = ['B', 'G', 'R', 'RE1', 'RE2', 'RE3', 'NIR', 'PAN']
    hammer_bands = ['PAN', 'HS0', 'HS1', 'HS2', 'HS3', 'HS4', 'HS5', 'HS6', 'HS7', 'HS8', 'HS9',
                    'HS10', 'HS11', 'HS12', 'HS13', 'HS14', 'HS15', 'HS16', 'HS17', 'HS18', 'HS19',
                    'HS20', 'HS21', 'HS22', 'HS23', 'HS24', 'HS25', 'HS26', 'HS27', 'HS28', 'HS29',
                    'HS30'] # 'TCI'

    if all([b in bands for b in menut_bands]):
        sensor = 'OpenCosmos_Menut'
        sensor_bands = [b for b in menut_bands]
    if all([b in bands for b in platero_bands]):
        sensor = 'OpenCosmos_Platero'
        sensor_bands = [b for b in platero_bands]
    if all([b in bands for b in hammer_bands]):
        sensor = 'OpenCosmos_Hammer'
        sensor_bands = [b for b in hammer_bands]
    if sensor is None: return

    ## list paths
    paths = {}
    paths['bands'] = {}
    for b in sensor_bands:
        for fi, f in enumerate(files):
            sp = basenames[fi].split('_')
            if len(sp) < 2: continue
            if basenames[fi].split('_')[-2] == b:
                paths['bands'][b] = f
                break

    ## text files
    text_files = ['ancillary', 'metadata', 'all',  'geoaccuracy_alignment_stats', 'matched', 'rpc']
    for k in text_files:
        file = [f for fi, f in enumerate(files) if basenames[fi].endswith('{}.json'.format(k))]
        if len(file) == 1: paths[k] = file[0]
        file = [f for fi, f in enumerate(files) if basenames[fi].endswith('{}.txt'.format(k))]
        if len(file) == 1: paths[k] = file[0]

    if ('ancillary' in paths) & ('metadata' in paths):
        return(sensor, paths)
