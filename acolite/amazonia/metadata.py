## def metadata
## read AMAZONIA1 metadata
##
## written by Quinten Vanhellemont, RBINS
## 2022-01-14
## modifications: 2022-02-06 (QV) added vza, attitude and ephemeris data
##                2022-02-16 (QV) fixed parsing of CBERS4A data which doesn't have camera tags

def metadata(metafile):
    import dateutil.parser
    from xml.dom import minidom

    import acolite as ac
    import numpy as np

    try:
        xmldoc = minidom.parse(metafile)
    except:
        print('Error opening metadata file.')
        sys.exit()

    metadata = {}
    for camera in ['leftCamera', 'rightCamera']:
        cmeta = {}
        ct = xmldoc.getElementsByTagName(camera)
        if len(ct) == 0:
            ##CBERS?
            camera = 'main'
            ct = xmldoc
        else:
            ## AMAZONIA
            ct = ct[0]
        if camera in cmeta: continue

        main_tag = 'orientationAngle'
        tags = ['degree', 'minute', 'second', 'millisecond']
        tags_out = ['oA_degree', 'oA_minute', 'oA_second', 'oA_millisecond']
        for t in ct.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    cmeta[tags_out[i]] = float(node[0].firstChild.nodeValue)
        cmeta['vza'] = cmeta['oA_degree'] + cmeta['oA_minute']/60 +  (cmeta['oA_second']+cmeta['oA_millisecond']/1000)/3600

        main_tag = 'satellite'
        tags = ['name', 'number', 'instrument']
        tags_out = ['name', 'number', 'instrument']
        for t in ct.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    cmeta[tags_out[i]] = node[0].firstChild.nodeValue

        main_tag = 'availableBands'
        for t in ct.getElementsByTagName(main_tag):
            nodes = t.getElementsByTagName('band')
            cmeta['bands'] = []
            for node in nodes:
                cmeta['bands'].append(node.firstChild.nodeValue)

        main_tag = 'image'
        tags = ['columns', 'lines']
        tags_out = ['columns', 'lines']
        for t in ct.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    cmeta[tags_out[i]] = int(node[0].firstChild.nodeValue)

        main_tag = 'sunPosition'
        tags = ["elevation", "sunAzimuth"]
        tags_out = ["sun_elevation", 'saa']
        for t in ct.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    cmeta[tags_out[i]] = float(node[0].firstChild.nodeValue)
        cmeta['sza'] = 90 - cmeta['sun_elevation']


        main_tag = 'viewing'
        tags = ['begin', 'center', 'end']
        tags_out = ['viewing_begin', 'viewing_center', 'viewing_end']
        for t in ct.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    cmeta[tags_out[i]] = node[0].firstChild.nodeValue


        main_tag = 'timeStamp'
        tags = ['begin', 'center', 'end']
        tags_out = ['timestamp_begin', 'timestamp_center', 'timestamp_end']
        for t in ct.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    cmeta[tags_out[i]] = node[0].firstChild.nodeValue


        band_tags = ['gain', 'integrationTime', 'absoluteCalibrationCoefficient']
        for main_tag in band_tags:
            for t in ct.getElementsByTagName(main_tag):
                nodes = t.getElementsByTagName('band')
                for node in nodes:
                    band = node.getAttribute('name')
                    cmeta['{}_{}'.format(band, main_tag)] = float(node.firstChild.nodeValue)


        ## find attitude data
        main_tag = 'attitude'
        tags = ['time', 'roll', 'pitch', 'yaw', 'deltaRoll', 'deltaPitch', 'deltaYaw']
        tags_out = ['time', 'roll', 'pitch', 'yaw', 'deltaRoll', 'deltaPitch', 'deltaYaw']
        tmp = {}
        for t in ct.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                if tag not in tmp: tmp[tag] = []
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    if tag == 'time':
                        tmp[tag].append(node[0].firstChild.nodeValue)
                    else:
                        tmp[tag].append(float(node[0].firstChild.nodeValue))
        cmeta['attitudes'] = tmp

        ## find ephemeris data
        main_tag = 'ephemeris'
        tags = ['time', 'x', 'y', 'z', 'vx', 'vy', 'vz']
        tags_out = ['time', 'x', 'y', 'z', 'vx', 'vy', 'vz']
        tmp = {}
        for t in ct.getElementsByTagName(main_tag):
            for i,tag in enumerate(tags):
                if tag not in tmp: tmp[tag] = []
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    if tag == 'time':
                        tmp[tag].append(node[0].firstChild.nodeValue)
                    else:
                        tmp[tag].append(float(node[0].firstChild.nodeValue))
        cmeta['ephemerides'] = tmp

        ## add current meta to camera metadata
        metadata[camera] = cmeta

    metadata['sensor'] = '{}{}_{}'.format(metadata[camera]['name'],
                                          metadata[camera]['number'],
                                          metadata[camera]['instrument'])

    metadata['isotime'] = metadata[camera]['timestamp_center']


    return(metadata)
