## def region_find
## finds geojson region file
## written by Quinten Vanhellemont, RBINS
## 2018-03-14
## modifications:
##                2021-03-02 (QV) added to acolite-gen, added data dir from acolite

def region_find(region, maxlev = 3, ext='geojson'):
    import glob
    import acolite as ac
    region_dir = ac.config['data_dir'] + '/Regions'

    lev = 0
    gf = []
    while (len(gf)==0) & (lev < maxlev):
        sstr = '{}/{}.{}'.format(region_dir, '*/'*lev+'{}'.format(region),ext)
        gf = glob.glob(sstr)
        lev+=1

    if len(gf) >= 1: gf=gf[0]

    return(gf)
