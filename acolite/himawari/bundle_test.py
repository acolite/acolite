## def bundle_test
## finds Himawari files
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-19
## modifications: 2025-05-22 (QV) use date+time to list bundles

def bundle_test(bundle):
    import os, glob

    if os.path.isdir(bundle):
        files = glob.glob('{}/HS_H0[8|9]_*.DAT'.format(bundle))
        files += glob.glob('{}/HS_H0[8|9]_*.DAT.bz2'.format(bundle))
        files.sort()
    else:
        ## if file is given only list files with matching date+time (bn[6:21])
        dn = os.path.dirname(bundle)
        bn = os.path.basename(bundle)
        files = glob.glob('{}/HS_H0[8|9]{}*{}'.format(dn, bn[6:21], os.path.splitext(bn)[1]))
        files.sort()

    ## test files
    fd = {}
    for file in files:
        bn = os.path.basename(file)
        bn, ex = bn[0:bn.find('.')], bn[bn.find('.'):]
        sp = bn.split('_')

        ## Full disk data
        #HS_Hnn_YYYYMMDD_hhmm_Bbb_FLDK_Rjj_Skkll.DAT
        if (len(sp) == 8) & (sp[0] == 'HS') & (sp[5] == 'FLDK'):
            platform = sp[1]
            date = sp[2]
            time = sp[3]
            band = sp[4]
            observation = sp[5]
            resolution = sp[6]
            segment = sp[7]

            ## use date time to sort files
            date_time = '{}_{}'.format(date, time)

            ## add parameters to dict
            if platform not in fd: fd[platform] = {}
            if date_time not in fd[platform]: fd[platform][date_time] = {}
            if band not in fd[platform][date_time]: fd[platform][date_time][band] = {}

            ## check segments
            if segment in fd[platform][date_time][band]:
                print('Segment {} already in {}'.format(segment, fd[platform][date_time][band][segment]))
            else:
                fd[platform][date_time][band][segment] = {'path': file, 'resolution': resolution, 'observation': observation}
        else:
            continue

    if len(fd) == 0: return
    return(fd)
