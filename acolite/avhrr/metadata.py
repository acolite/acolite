## def metadata
## parses XML metadata from AVHRR L1B/1C European Data Set bundle
## written by Quinten Vanhellemont, RBINS
## 2024-11-14
## modifications: 2024-11-16 (QV) get orbitdirection as well

def metadata(metafile):
    import os, sys
    from xml.dom import minidom

    if not os.path.isfile(metafile):
        print('Metadata file not found.')
        sys.exit()

    try:
        xmldoc = minidom.parse(metafile)
    except:
        print('Error opening metadata file.')
        sys.exit()

    metadata = {}
    ## get image information
    metadata_tags = ["eop:identifier", "eop:orbitDirection",
                     "gml:beginPosition", "gml:endPosition", "gml:timePosition",
                     ]

    for tag in metadata_tags:
        node = xmldoc.getElementsByTagName(tag)
        if len(node) > 0: metadata[tag] = node[0].firstChild.nodeValue

    ## tags
    tags = ["eop:Platform", "eop:Instrument"]
    sub_tags = ["eop:shortName", "eop:serialIdentifier",]
    for base_tag in tags:
        for t in xmldoc.getElementsByTagName(base_tag):
            for tag in sub_tags:
                    node = t.getElementsByTagName(tag)
                    if len(node) > 0: metadata['{}:{}'.format(base_tag, tag)] = node[0].firstChild.nodeValue

    ## convert
    metadata['sensor'] = '{}{}_{}'.format(metadata['eop:Platform:eop:shortName'],
                                          metadata['eop:Platform:eop:serialIdentifier'],
                                          metadata['eop:Instrument:eop:shortName'])
    metadata['start_date'] = metadata['gml:beginPosition']
    metadata['end_date'] = metadata['gml:endPosition']

    return(metadata)
