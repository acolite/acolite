## def safe_test
## lists files in given S2 safe directory and returns dict with band and file path
## written by Quinten Vanhellemont, RBINS
## 2017-04-18
## modifications:
##                  2018-04-17 (QV) added check for . files and continue on length of split failure
##                  2018-04-18 (QV) added check for . files in GRANULE
##                  2018-06-07 (QV) added check for jp2 for band files
##                  2024-06-04 (QV) added Sen2Cor bundle testing

def safe_test(file):
    import os

    files = os.listdir(file)
    datafiles = {}
    for i, fname in enumerate(files):
        tmp = fname.split('.')
        path = '{}/{}'.format(file,fname)

        ## scene metadata
        if (tmp[-1] == ('xml')) & ('MTD' in tmp[0]):
            datafiles['metadata'] = {"path":path, "fname":fname}

        ## granules
        if (fname == 'GRANULE'):
            granules = os.listdir(path)

            datafiles['granules'] = [] #granules
            for granule in granules:
                if granule[0]=='.':continue
                datafiles['granules'].append(granule)
                granule_data = {}
                split = granule.split('_')
                tile = split[1]
                date = split[-1]
                path = '{}/{}/{}/'.format(file,fname,granule)

                granule_files = os.listdir(path)
                for j, grfname in enumerate(granule_files):
                    tmp = grfname.split('.')
                    path = '{}/{}/{}/{}'.format(file,fname,granule,grfname)

                    ## scene metadata
                    if (tmp[-1] == ('xml')) & ('MTD' in tmp[0]):
                        granule_data['metadata'] = {"path":path,
                                                    "fname":grfname}
                    ## band files
                    if (grfname == 'IMG_DATA'):
                        if granule[0:3] == 'L2A':  ## sen2cor
                            res_dirs = os.listdir('{}/{}/{}/{}/'.format(file,fname,granule,grfname))
                            res_dirs.sort()
                            bands = []
                            for res_dir in res_dirs:
                                if res_dir in ['R10m', 'R20m', 'R60m']:
                                    bands += os.listdir('{}/{}/{}/{}/{}/'.format(file,fname,granule,grfname, res_dir))
                        else:
                            bands = os.listdir('{}/{}/{}/{}/'.format(file,fname,granule,grfname))

                        for band in bands:
                            if band[0] == '.': continue
                            if band[-3:] != 'jp2': continue
                            tmp=band.split('.')
                            ext = tmp[-1]
                            fsplit = tmp[0].split('_')
                            if len(fsplit) == 3: btile,btime,bid = fsplit
                            elif len(fsplit) == 11: btile,btime,bid = fsplit[-2],fsplit[-4],fsplit[-1]
                            elif len(fsplit) == 4: btile,btime,bid,bres = fsplit[0],fsplit[1],fsplit[2], fsplit[3] ## sen2cor
                            else:
                                print(fsplit)
                                continue
                            if granule[0:3] == 'L2A':  ## sen2cor
                                res_dir = 'R{}'.format(bres)
                                path = '{}/{}/{}/{}/{}/{}'.format(file,fname,granule,grfname,res_dir, band)
                            else:
                                path = '{}/{}/{}/{}/{}'.format(file,fname,granule,grfname,band)

                            band_name = 'B{}'.format(bid[1:3].lstrip('0'))
                            if band_name not in granule_data:
                                granule_data[band_name] = {"path":path,
                                                     "fname":band, 'tile':btile, 'time':btime, 'bid':bid, 'band_name':band_name}
                            else:
                                ## Sen2Cor provides bands at non native resolution as well
                                granule_data[band_name+'_'+bres] = {"path":path,
                                                     "fname":band, 'tile':btile, 'time':btime, 'bid':bid, 'band_name':band_name}

                    datafiles[granule] = granule_data
    return(datafiles)
