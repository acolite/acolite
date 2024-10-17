## def ocssw_l1a_to_l1b
## download L1A and make L1B and GEO files
## written by Quinten Vanhellemont, RBINS
## 2024-10-16
## modifications: 2024-10-16 (QV) adapted from viirs_l1a_to_l1b, currently HawkEye and MODIS
##                2024-10-17 (QV) added bzip2 extraction for MODIS

def ocssw_l1a_to_l1b(inputfile, local = None, override_download = False, override_conversion = False,
                     url_base = 'https://oceandata.sci.gsfc.nasa.gov/ob/getfile/'):
    import os, subprocess, time, bz2
    import acolite as ac
    try:
        ocsswroot = os.environ['OCSSWROOT']
    except:
        print('OCSSWROOT needs to be set as environment variable.')
        return

    wd = os.getcwd()
    bn = os.path.basename(inputfile)
    dn = os.path.dirname(inputfile)

    if ('.L1A.nc' not in bn) & ('L1A_LAC' not in bn) & ('L1A_GAC' not in bn):
        print('L1A to L1B processing of {} not supported.'.format(bn))
        return

    if ('HAWKEYE' in bn.upper()):
        stype = 'HawkEye'
    elif (bn.upper()[0] in ['A', 'T']):
        stype = 'MODIS'
    else:
        print('L1A to L1B processing of {} not supported.'.format(bn))
        return

    if local is None:
        if dn == '':
            local = wd
        else:
            local = dn

    if not os.path.exists(local): os.makedirs(local)
    os.chdir(local)

    local_scene = '{}/{}'.format(local, bn)
    print(local_scene)

    ## downloading
    if (override_download) | (not os.path.exists(local_scene)):
        print('Downloading {}'.format(bn))
        if stype == 'MODIS':
            local_scene_bz2 = local_scene+'.bz2'
            ac.shared.download_file(url_base+bn, local_scene_bz2, verbosity=5)
            print('Extracting {}'.format(os.path.basename(local_scene_bz2)))
            with open(local_scene, 'wb') as fo, bz2.BZ2File(local_scene_bz2, 'rb') as fi:
                for data in iter(lambda : fi.read(100 * 1024), b''): fo.write(data)
        else:
            ac.shared.download_file(url_base+bn, local_scene, verbosity=5)
    else:
        print('{} exists locally'.format(bn))

    ## find out processing and set up ocssw
    if stype == 'HawkEye':
        geo = '{}/{}'.format(local, bn.replace('.L1A.nc','.GEO.nc'))
        l1b = '{}/{}'.format(local, bn.replace('.L1A.nc','.L1B.nc'))
        proc_geo = '{} {} {}'.format('geolocate_hawkeye', local_scene, geo)
        proc_l1b = '{} ifile={} geofile={} ofile={}'.format('l1bgen_generic', local_scene, geo, l1b)
        ## do we need to run ocssw
        make_geo = (override_conversion) | (not os.path.exists(geo))
        make_l1b = (override_conversion) | (not os.path.exists(l1b))

    if stype == 'MODIS':
        geo = '{}/{}'.format(local, bn.replace('.L1A_LAC','.GEO'))
        l1b = '{}/{}'.format(local, bn.replace('.L1A_LAC','.L1B'))
        l1b1k = l1b+'_1KM'
        l1bhk = l1b+'_HKM'
        l1bqk = l1b+'_QKM'
        #proc_att = '{} --refreshDB {}'.format('modis_atteph', local_scene)
        proc_geo = '{} {} -o {}'.format('modis_GEO', local_scene, geo)
        proc_l1b = '{} -o {} -k {} -q {} {} {}'.format('modis_L1B', l1b1k, l1bhk, l1bqk, local_scene, geo)

        ## do we need to run ocssw
        make_geo = (override_conversion) | (not os.path.exists(geo))
        make_l1b = (override_conversion) | (not os.path.exists(l1b1k)) | \
                                           (not os.path.exists(l1bhk)) | \
                                           (not os.path.exists(l1bqk))

    if (make_geo) | (make_l1b):
        ## set up bash script
        lines = ['#!/bin/bash',
                 'export OCSSWROOT={};'.format(ocsswroot),
                 'source $OCSSWROOT/OCSSW_bash.env']
        ## add processing lines
        #if (make_geo) & (stype == 'MODIS'): lines.append(proc_att)
        if make_geo: lines.append(proc_geo)
        if make_l1b: lines.append(proc_l1b)
        with open('ocssw_run.sh', 'w') as f:
            for line in lines:
                f.write(line)
                #print(line)
                f.write('\n')
        ## run the bash script
        print('Running OCSSW processing')
        start = time.time()
        ret = subprocess.run('bash ocssw_run.sh', capture_output=True, shell=True)
        print("Finished OCSSW processing, elapsed Time: {:.1f}s".format(time.time() - start))
        #print(ret.stdout.decode())
    else:
        print('GEO and L1B files exist locally')

    os.chdir(wd)

    ret = {}
    if stype == 'HawkEye':
        ret['l1b'] = l1b
        ret['geo'] = geo
    if stype == 'MODIS':
        ret['l1b_1km'] = l1b1k
        ret['l1b_hkm'] = l1bhk
        ret['l1b_qkm'] = l1bqk
        ret['geo'] = geo

    if all(os.path.exists(ret[f]) for f in ret):
        return(ret)
    else:
        print('Error in OCSSW L1A to L1B conversion')
        print('Are the right OCSSW tools installed?')
        return
