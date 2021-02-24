## def init
## initialised interpolators for geolocation of Pleiades image
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-17
## modifications:
##

def init(metadata):
    from scipy.interpolate import interp2d

    ## Upper Left P1
    p1_lon=metadata['VERTICES']['UL']['LON']
    p1_lat=metadata['VERTICES']['UL']['LAT']
    p1_col=metadata['VERTICES']['UL']['COL']
    p1_row=metadata['VERTICES']['UL']['ROW']

    ## Upper Right P2
    p2_lon=metadata['VERTICES']['UR']['LON']
    p2_lat=metadata['VERTICES']['UR']['LAT']
    p2_col=metadata['VERTICES']['UR']['COL']
    p2_row=metadata['VERTICES']['UR']['ROW']

    ## Lower Right P3
    p3_lon=metadata['VERTICES']['LR']['LON']
    p3_lat=metadata['VERTICES']['LR']['LAT']
    p3_col=metadata['VERTICES']['LR']['COL']
    p3_row=metadata['VERTICES']['LR']['ROW']

    ## Lower Left P4
    p4_lon=metadata['VERTICES']['LL']['LON']
    p4_lat=metadata['VERTICES']['LL']['LAT']
    p4_col=metadata['VERTICES']['LL']['COL']
    p4_row=metadata['VERTICES']['LL']['ROW']

    ## Center PC
    pc_lon=metadata['VERTICES']['C']['LON']
    pc_lat=metadata['VERTICES']['C']['LAT']
    pc_col=metadata['VERTICES']['C']['COL']
    pc_row=metadata['VERTICES']['C']['ROW']

    ncols = float(metadata['NCOLS'])
    nrows = float(metadata['NROWS'])

    if True:
        ## get vertex image location
        pcol = [p1_col, p2_col, p3_col, p4_col]
        prow = [p1_row, p2_row, p3_row, p4_row]

        ## get vertex geolocation
        plon = [p1_lon, p2_lon, p3_lon, p4_lon]
        plat = [p1_lat, p2_lat, p3_lat, p4_lat]
    else:
        ## get vertex image location
        pcol = [p1_col, p2_col, p3_col, p4_col, pc_col]
        prow = [p1_row, p2_row, p3_row, p4_row, pc_row]

        ## get vertex geolocation
        plon = [p1_lon, p2_lon, p3_lon, p4_lon, pc_lon]
        plat = [p1_lat, p2_lat, p3_lat, p4_lat, pc_lat]

    ## set up interpolator
    zlon = interp2d(pcol, prow, plon)
    zlat = interp2d(pcol, prow, plat)

    ## set up interpolator for lat/lon -> row/col
    zcol = interp2d(plon, plat, pcol)
    zrow = interp2d(plon, plat, prow)

    return (zlon, zlat), (zcol, zrow)
