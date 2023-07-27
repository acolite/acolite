## def viirs_l1a_to_l1b
## download VIIRS L1A and make L1B and GEO files
## written by Quinten Vanhellemont, RBINS
## 2023-04-01
## modifications: 2023-04-02 (QV) functionised

def viirs_l1a_to_l1b(scene, local = None, override_download = False, override_conversion = False,
                     url_base = 'https://oceandata.sci.gsfc.nasa.gov/ob/getfile/'):
    import os, subprocess, time
    import acolite as ac
    try:
        ocsswroot = os.environ['OCSSWROOT']
    except:
        print('OCSSWROOT needs to be set as environment variable.')
        return

    wd = os.getcwd()
    bn = os.path.basename(scene)
    dn = os.path.dirname(scene)

    if ('VIIRS' not in scene) | ('.L1A.nc' not in scene):
        print('VIIRS L1A to L1B processing of {} not supported.'.format(scene))
        return

    if local is None: local = wd
    if not os.path.exists(local): os.makedirs(local)
    os.chdir(local)

    local_scene = '{}/{}'.format(local, scene)

    ## downloading
    if (override_download) | (not os.path.exists(local_scene)):
        print('Downloading {}'.format(scene))
        ac.shared.download_file(url_base+scene, local_scene, verbosity=5)
    else:
        print('{} exists locally'.format(scene))

    ## find out processing
    l1b_mod = scene.replace('.L1A.nc','.L1B_MOD.nc')
    l1b_img = scene.replace('.L1A.nc','.L1B_IMG.nc')
    geo_mod = scene.replace('.L1A.nc','.GEO_MOD.nc')
    geo_img = scene.replace('.L1A.nc','.GEO_IMG.nc')

    ## set up ocssw code
    proc_geo = '{} ifile={} geofile_img={} geofile_mod={}'.format('geolocate_viirs', scene,
                                                                           geo_img, geo_mod)
    proc_l1b = '{} ifile={} l1bfile_img={} l1bfile_mod={}'.format('calibrate_viirs', scene,
                                                                           l1b_img, l1b_mod)
    ## do we need to run ocssw
    make_geo = (override_conversion) | (not os.path.exists(geo_mod)) | (not os.path.exists(geo_img))
    make_l1b = (override_conversion) | (not os.path.exists(l1b_mod)) | (not os.path.exists(l1b_img))
    if (make_geo) | (make_l1b):
        ## set up bash script
        lines = ['#!/bin/bash',
                 'export OCSSWROOT={};'.format(ocsswroot),
                 'source $OCSSWROOT/OCSSW_bash.env']
        ## add processing lines
        if make_geo: lines.append(proc_geo)
        if make_l1b: lines.append(proc_l1b)
        with open('viirs_run.sh', 'w') as f:
            for line in lines:
                f.write(line)
                #print(line)
                f.write('\n')
        ## run the bash script
        print('Running ocssw VIIRS processing')
        start = time.time()
        ret = subprocess.run('bash viirs_run.sh', capture_output=True, shell=True)
        print("Finished ocssw VIIRS processing, elapsed Time: {:.1f}s".format(time.time() - start))
        #print(ret.stdout.decode())
    else:
        print('GEO and L1B files exist locally')

    os.chdir(wd)

    ret = {}
    ret['l1b_mod'] = '{}/{}'.format(local,l1b_mod)
    ret['geo_mod'] = '{}/{}'.format(local,geo_mod)
    ret['l1b_img'] = '{}/{}'.format(local,l1b_img)
    ret['geo_img'] = '{}/{}'.format(local,geo_img)
    if all(os.path.exists(ret[f]) for f in ret):
        return(ret)
    else:
        print('Error in VIIRS L1A to L1B conversion')
        print('Are OCSSW viirs_l1_bin and viirs_l1 tools installed?')
        return
