## def metadata_parse
## parses XML metadata from Huanjing bundle
## written by Quinten Vanhellemont, RBINS
## 2025-07-29
## modifications:

def metadata_parse(metafile):
    import os, sys, fnmatch, dateutil.parser
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
    metadata_tags = [\
    'SatelliteID', 'ReceiveStationID', 'SensorID', 'ReceiveTime',
    'OrbitID', 'OrbitType','AttType', 'StripID',
    'ProduceType', 'SceneID', 'DDSFlag', 'ProductID',
    'ProductLevel', 'ProductFormat','ProduceTime', 'Bands',
    'ScenePath','SceneRow','SatPath','SatRow', 'SceneCount','SceneShift',
    'StartTime','EndTime',
    'CenterTime', 'StartLine', 'EndLine',
    'ImageGSD','WidthInPixels','HeightInPixels','WidthInMeters','HeightInMeters',
    'RegionName','CloudPercent','DataSize',
    'RollViewingAngle','PitchViewingAngle',
    'PitchSatelliteAngle','RollSatelliteAngle','YawSatelliteAngle',
    'SolarAzimuth', 'SolarZenith', 'SatelliteAzimuth','SatelliteZenith',
    'GainMode', 'IntegrationTime', 'IntegrationLevel',
    'MapProjection', 'EarthEllipsoid', 'ZoneNo',
    'ResamplingKernel', 'HeightMode', 'EphemerisData', 'AttitudeData', 'RadiometricMethod',
    'MtfCorrection','Denoise','RayleighCorrection','UsedGCPNo'
    'CenterLatitude', 'CenterLongitude',
    'TopLeftLatitude', 'TopLeftLongitude', 'TopRightLatitude', 'TopRightLongitude',
    'BottomRightLatitude', 'BottomRightLongitude', 'BottomLeftLatitude', 'BottomLeftLongitude',
    'TopLeftMapX','TopLeftMapY','TopRightMapX','TopRightMapY',
    'BottomRightMapX','BottomRightMapY','BottomLeftMapX','BottomLeftMapY',
    'DataFile','BrowseFile','ThumbFile']

    for tag in metadata_tags:
        node = xmldoc.getElementsByTagName(tag)
        if len(node) > 0:
            if node[0].firstChild is None: continue
            metadata[tag] = node[0].firstChild.nodeValue

    return(metadata)
