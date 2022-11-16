## def metadata_parse
## parse EnMAP metadata
## written by Quinten Vanhellemont, RBINS
## 2022-09-19
## modifications:

def metadata_parse(metafile, parse_state_vectors=False):
    import dateutil.parser
    from xml.dom import minidom
    import numpy as np

    try:
        xmldoc = minidom.parse(metafile)
    except:
        print('Error opening metadata file.')
        sys.exit()

    xml_main = xmldoc.firstChild

    metadata = {}

    ## required tags
    tags = ['sunAzimuthAngle', 'sunElevationAngle',
            'acrossOffNadirAngle', 'alongOffNadirAngle','sceneAzimuthAngle',
            'mapProjection', 'startTime', 'stopTime', 'datatakeStart', 'datatakeStop',
            'heightOfOrthoScene', 'widthOfOrthoScene',
            'sceneSZA', 'sceneAOT', 'sceneWV', 'projection',
            'format', 'level', 'mission', 'sensor',
            ]

    for tag in tags:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) == 0: continue

        if tdom[0].hasChildNodes():
            nodes = tdom[0].childNodes
            cc = tdom[0].firstChild
            while cc is not None:
                if cc.firstChild is not None:
                    name = cc.nodeName
                    val = cc.firstChild.nodeValue.strip()
                else:
                    name = None
                    val = cc.nodeValue.strip()
                if val not in ['', '\n']:
                    try:
                        val = float(val)
                    except:
                        pass
                    if name is None:
                        metadata[tag] = val
                    else:
                        if tag not in metadata: metadata[tag] = {}
                        metadata[tag][name] = val
                cc = cc.nextSibling

    ## parse band information
    band_data = {}
    tag = 'bandID'
    tdom = xmldoc.getElementsByTagName(tag)
    for t in tdom:
        band_number  = t.getAttribute('number')
        if band_number not in band_data: band_data[band_number] = {}
        if t.hasChildNodes():
            nodes = t.childNodes
            cc = t.firstChild
            while cc is not None:
                if cc.firstChild is not None:
                    name = cc.nodeName
                    val = cc.firstChild.nodeValue.strip()
                    try:
                        val = float(val)
                    except:
                        pass
                    band_data[band_number][name] = val
                cc = cc.nextSibling

    ## parse state vector information
    if parse_state_vectors:
        state_vectors = {}
        tag = 'stateVector'
        tdom = xmldoc.getElementsByTagName(tag)
        for t in tdom:
            svnum  = t.getAttribute('num')
            svmaneuver = t.getAttribute('maneuver')
            svqualInd = t.getAttribute('qualInd')
            if svnum not in state_vectors:
                state_vectors[svnum] = {'num': svnum, 'maneuver': svmaneuver, 'qualInd': svqualInd}
            if t.hasChildNodes():
                nodes = t.childNodes
                cc = t.firstChild
                while cc is not None:
                    if cc.firstChild is not None:
                        name = cc.nodeName
                        val = cc.firstChild.nodeValue.strip()
                        try:
                            val = float(val)
                        except:
                            pass
                        state_vectors[svnum][name] = val
                    cc = cc.nextSibling

        return(metadata, band_data, state_vectors)
    else:
        return(metadata, band_data)
