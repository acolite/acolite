## def metadata
## reads some metadata from DIMAP dim file
## written by Quinten Vanhellemont, RBINS
## 2023-02-14

def metadata(dimfile):
    import xml.etree.ElementTree as ET

    ## set up meta dict
    meta = {}

    ## open dim file
    tree = ET.parse(dimfile)
    #all_elem = list(tree.iter())
    prod = tree.find('Production')
    meta['start_date'] = prod.find('PRODUCT_SCENE_RASTER_START_TIME').text
    meta['end_date'] = prod.find('PRODUCT_SCENE_RASTER_STOP_TIME').text

    mt = tree.find('Metadata_Id')
    meta['met_format'] = mt.find('METADATA_FORMAT').text

    mt = tree.find('Dataset_Id')
    meta['met_dataset'] = mt.find('DATASET_SERIES').text

    ## parse TPG info
    tpg = tree.find('Tie_Point_Grids')
    tpg_atts = {}
    if tpg is not None:
        for t in tpg.findall('Tie_Point_Grid_Info'):
            tpg_att = {k.tag : k.text for k in t.iter()}
            tpg_atts[tpg_att['TIE_POINT_GRID_NAME']] = tpg_att
    tree = None
    meta['tpg_atts'] = tpg_atts
    return(meta)
