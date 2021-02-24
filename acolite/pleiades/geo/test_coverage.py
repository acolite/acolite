## def test_coverage
## tests if ROI given by limit [S,W,N,E] is covered by given image
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-17
## modifications: QV 2017-03-23 added quick and dirty subset to returned file lists
##                2020-10-18 (QV) adapted for acg

def test_coverage(d, limit, verbose=False):
    import acolite as ac

    if type(d) == dict:
        meta = d
    else:
        ## check if we are dealing with Pléiades image bundle
        ifile,mfile = ac.pleiades.bundle_test(d, listpan=False)
        ifile = ifile[0]
        mfile = mfile[0]
        ## parse metadata
        meta = ac.pleiades.metadata_parse(mfile)

    ## check for coverage
    sub = ac.pleiades.geo.crop(meta, limit)

    ncols = float(meta['NCOLS'])
    nrows = float(meta['NROWS'])

    yu=0
    yl=0
    xl=0
    xr=0

    ### check longitudes
    if sub[0] >= ncols:
        if verbose: print('Western longitude out of eastern bound of scene.')
        xl=1

    if (sub[0]+sub[2]) >= ncols:
        if verbose: print('Eastern longitude out of eastern bound of scene.')
        xr=1

    if sub[0] <= 1:
        if verbose: print('Eastern longitude out of western bound of scene.')
        xl=-1

    if (sub[0]+sub[2]) <= 1:
        if verbose: print('Western longitude out of western bound of scene.')
        xr=-1

    ## check latitudes
    if sub[1] <= 1:
        if verbose: print('Northern longitude out of northern bound of scene.')
        yu=1

    if (sub[1]+sub[3]) <= 1:
        if verbose: print('Southern longitude out of northern bound of scene.')
        yl=1

    if sub[1] >= nrows:
        if verbose: print('Northern longitude out of southern bound of scene.')
        yu=-1

    if (sub[1]+sub[3]) >= nrows:
        if verbose: print('Southern longitude out of southern bound of scene.')
        yl=-1

    ## check dimensions
    if (sub[2] == 1) or (sub[3] == 1):
        if verbose: print('Crop dimensions of 1: region probably out of scene.')
        yl=-1
        yu=-1
        xl=-1
        xr=-1

    if (yl == 0) & (yu == 0) & (xl == 0) & (xr == 0):
        if verbose: print('Region fully in scene.')
        return(False)

    if (xl == 1):
        if verbose: print('Region fully west of scene')
        return(True)

    if (xr == -1):
        if verbose: print('Region fully east of scene')
        return(True)

    if (yl == 1):
        if verbose: print('Region fully north of scene')
        return(True)

    if (yu == -1):
        if verbose: print('Region fully south of scene')
        return(True)

    if verbose: print('Region partially within bounds of scene')
    return(False)
