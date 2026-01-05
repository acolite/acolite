## def bundle_test
## identifies files for Hyperfield data
## written by Quinten Vanhellemont, RBINS
## 2026-01-06
## modifications:

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
        print('Error')
        return

    basenames = [os.path.basename(f) for f in files]

    sensor = None
    level = None
    paths = {}
    for file in files:
        f0 = os.path.basename(os.path.splitext(file)[0])
        if file.endswith('.json'):
            sp = f0.split('_')
            if sp[1] == 'L1C':
                level = sp[1]
                sensor = sp[0]
                paths['metadata'] = file

        if file.endswith('.tif'):
            paths[f0] = file

    if sensor is None: return
    if ('L1C' in paths) & ('metadata' in paths):
        return(sensor, paths)
