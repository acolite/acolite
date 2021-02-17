## def metadata_granule
## imports S2 granule metadata
## written by Quinten Vanhellemont, RBINS
## 2017-04-18
## modifications:
##                2018-07-18 (QV) changed acolite import name
##                2020-10-28 (QV) fill nans in angles grids
##                2021-02-11 (QV) adapted for acolite-gen, renamed from granule_meta
##                2021-02-17 (QV) added fillnan keyword, added per detector grids, renamed safe_tile_grid

def metadata_granule(metafile, fillnan=False):
    import dateutil.parser
    from xml.dom import minidom

    import acolite as ac
    import numpy as np
    import copy

    try:
        xmldoc = minidom.parse(metafile)
    except:
        print('Error opening metadata file.')
        sys.exit()

    xml_main = xmldoc.firstChild

    metadata = {}
    tags = ['TILE_ID','DATASTRIP_ID', 'SENSING_TIME']
    for tag in tags:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) > 0: metadata[tag] = tdom[0].firstChild.nodeValue

    #Geometric_Info
    grids = {'10':{}, '20':{}, '60':{}}
    Geometric_Info = xmldoc.getElementsByTagName('n1:Geometric_Info')
    if len(Geometric_Info) > 0:
        for tg in Geometric_Info[0].getElementsByTagName('Tile_Geocoding'):
            tags = ['HORIZONTAL_CS_NAME','HORIZONTAL_CS_CODE']
            for tag in tags:
                tdom = tg.getElementsByTagName(tag)
                if len(tdom) > 0: metadata[tag] = tdom[0].firstChild.nodeValue

            for sub in  tg.getElementsByTagName('Size'):
                res = sub.getAttribute('resolution')
                grids[res]['RESOLUTION'] = float(res)
                tags = ['NROWS','NCOLS']
                for tag in tags:
                    tdom = sub.getElementsByTagName(tag)
                    if len(tdom) > 0: grids[res][tag] = int(tdom[0].firstChild.nodeValue)

            for sub in  tg.getElementsByTagName('Geoposition'):
                res = sub.getAttribute('resolution')
                tags = ['ULX','ULY','XDIM','YDIM']
                for tag in tags:
                    tdom = sub.getElementsByTagName(tag)
                    if len(tdom) > 0: grids[res][tag] = int(tdom[0].firstChild.nodeValue)

        for ta in Geometric_Info[0].getElementsByTagName('Tile_Angles'):
            ## sun angles
            sun_angles={}
            for tag in ['Zenith','Azimuth']:
                for sub in ta.getElementsByTagName('Sun_Angles_Grid')[0].getElementsByTagName(tag):
                    sun_angles[tag] = ac.sentinel2.grid_geom(sub)

            for sub in ta.getElementsByTagName('Mean_Sun_Angle'):
                sun_angles['Mean_Zenith'] = float(sub.getElementsByTagName('ZENITH_ANGLE')[0].firstChild.nodeValue)
                sun_angles['Mean_Azimuth'] = float(sub.getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.nodeValue)

            ## view angles
            view_angles={} # merged detectors
            view_angles_det={} # separate detectors
            for sub in ta.getElementsByTagName('Viewing_Incidence_Angles_Grids'):
                band = sub.getAttribute('bandId')
                detector = sub.getAttribute('detectorId')
                if band not in view_angles_det: view_angles_det[band]={}
                if detector not in view_angles_det[band]: view_angles_det[band][detector]={}
                band_view = {}
                for tag in ['Zenith','Azimuth']:
                    ret = ac.sentinel2.grid_geom(sub.getElementsByTagName(tag)[0])
                    band_view[tag] = copy.copy(ret)
                    view_angles_det[band][detector][tag] = copy.copy(ret)
                if band not in view_angles.keys():
                    view_angles[band] = band_view
                else:
                    for tag in ['Zenith','Azimuth']:
                        mask = np.isfinite(band_view[tag]) & np.isnan(view_angles[band][tag])
                        view_angles[band][tag][mask] = band_view[tag][mask]

            if fillnan:
                for b,band in enumerate(view_angles.keys()):
                    for tag in ['Zenith','Azimuth']:
                        view_angles[band][tag] = ac.shared.fillnan(view_angles[band][tag])

            ## average view angle grid
            ave = {}
            for b,band in enumerate(view_angles.keys()):
                for tag in ['Zenith','Azimuth']:
                    data = view_angles[band][tag]
                    count = np.isfinite(data)*1
                    if b == 0:
                        ave[tag] = data
                        ave['{}_Count'.format(tag)] = count
                    else:
                        ave[tag] += data
                        ave['{}_Count'.format(tag)] += count
            for tag in ['Zenith','Azimuth']: view_angles['Average_View_{}'.format(tag)] = ave[tag] / ave['{}_Count'.format(tag)]

    metadata["GRIDS"] = grids
    metadata["VIEW"] = view_angles
    metadata["VIEW_DET"] = view_angles_det
    metadata["SUN"] = sun_angles

    ## some interpretation
    #metadata['TIME'] = dateutil.parser.parse(metadata['SENSING_TIME'])
    #metadata["DOY"] = metadata["TIME"].strftime('%j')
    #metadata["SE_DISTANCE"] = ac.shared.distance_se(metadata['DOY'])

    return(metadata)
