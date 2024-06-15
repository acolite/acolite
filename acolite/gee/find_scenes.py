## def find scenes
## finds L1 Landsat or Sentinel-2 scenes on GEE for a given limit ROI or lat/lon points
## written by Quinten Vanhellemont, RBINS
## 2022-04-12
## modifications: 2022-12-25 (QV) added SR option
##                2023-01-05 (QV) added L4 and L1-5 MSS
##                2024-06-15 (QV) change to S2_HARMONIZED

def find_scenes(isodate_start, isodate_end=None, day_range=1,
                surface_reflectance=False,
                limit=None, st_lat=None, st_lon=None, filter_tiles=None,
                sensors=['L4_TM', 'L5_TM', 'L7_ETM', 'L8_OLI', 'L9_OLI', 'S2A_MSI', 'S2B_MSI']):
    import ee
    #ee.Authenticate() ## assume ee use is authenticated in current environment
    #ee.Initialize()

    import dateutil.parser, datetime

    if filter_tiles is not None:
        if type(filter_tiles) is not list:
            filter_tiles = [filter_tiles]

    ## check isodate
    if isodate_start is None:
        print('Please provide start date.')
        return()
    else:
        dstart = dateutil.parser.parse(isodate_start)
        isodate_start = dstart.isoformat()[0:10]

    ## get date range
    if isodate_start == isodate_end: isodate_end = None
    if isodate_end is None:
        dend = dstart + datetime.timedelta(days=0)
    else:
        if isodate_end in ['now', 'today']:
            dend = datetime.datetime.now()
        else:
            dend = dateutil.parser.parse(isodate_end)
    dend += datetime.timedelta(days=1) ## add one day so end date is included
    isodate_end = dend.isoformat()[0:10]

    print('Date range {} {}'.format(isodate_start, isodate_end))

    ## identify collections
    collections = []
    landsats = []
    ## MultiSpectral Scanners
    if 'L1_MSS' in sensors: landsats.append('LM01')
    if 'L2_MSS' in sensors: landsats.append('LM02')
    if 'L3_MSS' in sensors: landsats.append('LM03')
    if 'L4_MSS' in sensors: landsats.append('LM04')
    if 'L5_MSS' in sensors: landsats.append('LM05')

    ## newer sensors
    if 'L4_TM' in sensors: landsats.append('LT04')
    if 'L5_TM' in sensors: landsats.append('LT05')
    if 'L7_ETM' in sensors: landsats.append('LE07')
    if 'L8_OLI' in sensors: landsats.append('LC08')
    if 'L9_OLI' in sensors: landsats.append('LC09')
    landsat_tiers = ['T1', 'T2']
    landsat_collections = ['C02']

    for landsat in landsats:
        for tier in landsat_tiers:
            for coll in landsat_collections:
                if surface_reflectance:
                    if landsat[1] == 'M':
                        print('No SR for MSS.')
                    else:
                        collections.append('{}/{}/{}/{}_L2'.format('LANDSAT', landsat, coll, tier))
                else:
                    if landsat[1] == 'M':
                        collections.append('{}/{}/{}/{}'.format('LANDSAT', landsat, coll, tier))
                    else:
                        collections.append('{}/{}/{}/{}_TOA'.format('LANDSAT', landsat, coll, tier))

    if ('S2A_MSI' in sensors) or ('S2B_MSI' in sensors):
        ## harmonized has scenes from new processing shifted to old processing
        ## we take the offset into account in agh for >= PB4 data
        if surface_reflectance:
            #collections += ['COPERNICUS/S2_SR'] # COPERNICUS/S2_SR_HARMONIZED
            collections += ['COPERNICUS/S2_SR_HARMONIZED'] # COPERNICUS/S2_SR superseded by COPERNICUS/S2_SR_HARMONIZED in Jun 2024
        else:
            #collections.append('COPERNICUS/S2') # 'COPERNICUS/S2_HARMONIZED'
            collections.append('COPERNICUS/S2_HARMONIZED') # COPERNICUS/S2 superseded by COPERNICUS/S2_HARMONIZED in Jun 2024

    print('Checking collections {}'.format(' '.join(collections)))

    ## set up region
    if limit is not None:
        region = ee.Geometry.BBox(limit[1], limit[0], limit[3], limit[2])
    elif (st_lon is not None) & (st_lat is not None):
        region = ee.Geometry.Point([st_lon, st_lat])
    else:
        print('Warning! No limit or st_lat, st_lon combination specified. Function may return too many images.')
        region = None

    ## set up ee date
    sdate=ee.Date(isodate_start)
    edate=ee.Date(isodate_end)

    ## search ee collections
    imColl = None
    for coll in collections:
        imC = ee.ImageCollection(coll).filterDate(sdate, edate)
        if region is not None: imC = imC.filterBounds(region)

        if imColl is None:
            imColl = imC
        else:
            imColl = imColl.merge(imC)

    iml = imColl.getInfo()
    nimages = len(iml['features'])

    images = []
    if nimages > 0:
        limages = imColl.toList(nimages).getInfo()
        for im in limages:
            if 'PRODUCT_ID' in im['properties']: ## Sentinel-2 image
                fkey = 'PRODUCT_ID'
                pid = im['properties'][fkey]
            elif 'LANDSAT_PRODUCT_ID' in im['properties']: ## Landsat image
                fkey = 'LANDSAT_PRODUCT_ID'
                pid = im['properties'][fkey]
            else: continue

            skip = False
            if filter_tiles is not None:
                skip = True
                for tile in filter_tiles:
                    if tile in pid: skip = False
            if skip: continue
            images.append((fkey,pid))
    return(images, imColl)
