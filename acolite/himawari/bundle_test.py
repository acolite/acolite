## def bundle_test
## finds Himawari files
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-19
## modifications:

def bundle_test(bundle):
    import os, glob

    if os.path.isdir(bundle):
        files = glob.glob('{}/HS_H0[8|9]_*.DAT'.format(bundle))
        files += glob.glob('{}/HS_H0[8|9]_*.DAT.bz2'.format(bundle))
        files.sort()
    else:
        dn = os.path.dirname(bundle)
        bn = os.path.basename(bundle)
        files = glob.glob('{}/HS_H0[8|9]{}*{}'.format(dn, bn[6:16], os.path.splitext(bn)[1]))
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

            ## add parameters to dict
            if platform not in fd: fd[platform] = {}
            if date not in fd[platform]: fd[platform][date] = {}
            if band not in fd[platform][date]: fd[platform][date][band] = {}

            ## check segments
            if segment in fd[platform][date][band]:
                print('Segment {} already in {}'.format(segment, fd[platform][date][band][segment]))
            else:
                fd[platform][date][band][segment] = {'path': file, 'resolution': resolution, 'observation': observation}
        else:
            continue

    if len(fd) == 0: return
    return(fd)
