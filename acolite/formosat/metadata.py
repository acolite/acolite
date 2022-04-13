## def metadata
## read FORMOSAT5 metadata
##
## written by Quinten Vanhellemont, RBINS
## 2022-04-12
## modifications:

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

    meta = {}
    tags = ['DATASET_NAME', 'MISSION', 'MISSION_INDEX', 'INSTRUMENT']
    tags += ['IMAGING_DATE', 'IMAGING_TIME', 'SUN_AZIMUTH', 'SUN_ELEVATION',
             'VIEWING_ANGLE', 'VIEWING_ANGLE_CROSS_TRACK','VIEWING_ANGLE_ALONG_TRACK',
             'SATELLITE_INCIDENCE_ANGLE','SATELLITE_AZIMUTH_ANGLE']

    for tag in tags:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) > 0:
            meta[tag] = tdom[0].firstChild.nodeValue

    tag = 'DATA_FILE_PATH'
    tdom = xmldoc.getElementsByTagName(tag)
    if len(tdom) > 0:
        meta[tag] = tdom[0].getAttribute('href')


    meta['isodate']='{}T{}'.format(meta['IMAGING_DATE'],meta['IMAGING_TIME'])
    meta['saa']=float(meta['SUN_AZIMUTH'])
    meta['sza']=90-float(meta['SUN_ELEVATION'])
    meta['vza']=float(meta['VIEWING_ANGLE'])
    meta['vaa']=float(meta['SATELLITE_AZIMUTH_ANGLE'])

    meta['satellite'] = '{}{}'.format(meta['MISSION'], meta['MISSION_INDEX'])
    meta['sensor']='{}_{}'.format(meta['satellite'], meta['INSTRUMENT'])



    tdom = xmldoc.getElementsByTagName('Spectral_Band_Info')
    if len(tdom) > 0:
        tdom[0]

    band_tags = ['BAND_INDEX', 'GAIN_NUMBER', 'PHYSICAL_CALIBRATION_DATE', 'PHYSICAL_UNIT',
                'PHYSICAL_GAIN', 'PHYSICAL_BIAS']
    band_names = ['Blue', 'Green', 'Red', 'NIR']

    bands = {}
    for t in tdom:
        band_dct = {}
        for tag in band_tags:
            bt = t.getElementsByTagName(tag)
            if len(bt) > 0:
                band_dct[tag]=bt[0].firstChild.nodeValue
        bands[band_dct['BAND_INDEX']] = band_dct

    for bi, b in enumerate(band_names):
        bands[str(bi+1)]['NAME'] = b

    for b in bands:
        for k in bands[b]:
            if k in ['PHYSICAL_GAIN', 'PHYSICAL_BIAS']:
                bands[b][k] = float(bands[b][k])
            meta['{}_{}'.format(bands[b]['NAME'], k)] = bands[b][k]

    return(meta)
