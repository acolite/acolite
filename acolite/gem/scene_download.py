## def scene_download
## download scene from GCS using gsutil
## written by Quinten Vanhellemont, RBINS
## 2021-03-11
## modifications:

def scene_download(im, output_base=None, download=True, verbosity=0):

    import os, subprocess
    if output_base is None: output_base = os.getcwd()

    if 'PRODUCT_ID' in im['properties']:
        pid = im['properties']['PRODUCT_ID']
        tile = pid.split('_')[-2]
        source = "gs://gcp-public-data-sentinel-2/tiles/{}/{}/{}/{}.SAFE".format(tile[1:3],tile[3],tile[4:6], pid)
        target = "'{}/S2/{}.SAFE'".format(output_base, pid)
    if 'LANDSAT_PRODUCT_ID' in im['properties']:
        pid = im['properties']['LANDSAT_PRODUCT_ID']
        if im['properties']['SPACECRAFT_ID'] == 'LANDSAT_5': sengcs = 'LT05'
        if im['properties']['SPACECRAFT_ID'] == 'LANDSAT_7': sengcs = 'LE07'
        if im['properties']['SPACECRAFT_ID'] == 'LANDSAT_8': sengcs = 'LC08'
        sensor = 'L{}'.format(im['properties']['SPACECRAFT_ID'][8:])
        coll = '01'
        path = '{}'.format(im['properties']['WRS_PATH']).zfill(3)
        row = '{}'.format(im['properties']['WRS_ROW']).zfill(3)
        source = "gs://gcp-public-data-landsat/{}/{}/{}/{}/{}/".format(sengcs, coll, path, row, pid)
        target = '{}/{}/{}/'.format(output_base,sensor,pid)

    #if not os.path.exists(target): os.makedirs(target)
    #cmd = "gsutil -m rsync -r {} {}".format(source, target)
    cmd = "gsutil rsync -r {} '{}'".format(source, target)

    if verbosity > 1: print(cmd)
    if download:
        if verbosity > 0: print('Syncing {} to {}'.format(pid, output_base))
        p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    return(target)
