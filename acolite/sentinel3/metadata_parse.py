## def metadata_parse
## parse sensor from Sentinel-3/OLCI or ENVISAT/MERIS xfdumanifest.xml file
## written by Quinten Vanhellemont, RBINS
## 2022-11-03
## modifications:
##

def metadata_parse(metafile):
    from xml.dom import minidom
    try:
        xmldoc = minidom.parse(metafile)
    except:
        print('Error opening metadata file.')

    metadata = {}
    tags = ['sentinel-safe:familyName', 'sentinel-safe:number',
            'envisat:productName', 'envisat:productType',
            'sentinel3:productName', 'sentinel3:productType']
    for tag in tags:
        node = xmldoc.getElementsByTagName(tag)
        for n in node:
            ptag = n.parentNode.tagName
            if ptag in ['sentinel-safe:platform', 'sentinel-safe:instrument']:
                if ptag not in metadata: metadata[ptag] = {}
                metadata[ptag][tag] = n.firstChild.nodeValue
            else:
                metadata[tag] = n.firstChild.nodeValue

    sensor = None
    if (metadata['sentinel-safe:platform']['sentinel-safe:familyName'] == 'Sentinel-3') &\
       (metadata['sentinel-safe:instrument']['sentinel-safe:familyName'] == 'Ocean Land Colour Instrument'):
        sensor = 'S3{}_OLCI'.format(metadata['sentinel-safe:platform']['sentinel-safe:number'])
    if (metadata['sentinel-safe:platform']['sentinel-safe:familyName'] == 'Envisat') &\
       (metadata['sentinel-safe:instrument']['sentinel-safe:familyName'] == 'MEdium Resolution Imaging Spectrometer'):
        sensor = 'EN1_MERIS'
    metadata['sensor'] = sensor

    return(metadata)
