## def polygon_from_wkt
## creates temporary json file if wkt polygon is given
## written by Quinten Vanhellemont, RBINS
## 2023-02-06
## modifications: 2023-02-14 (QV) return None if failed
##                2024-04-08 (QV) added geojson
##                2025-01-30 (QV) added encoding to output json file

def polygon_from_wkt(wkt, file=None):
    import os, json
    from osgeo import ogr
    import acolite as ac

    if file == None: file = ac.config['scratch_dir'] + '/polygon.json'
    odir = os.path.dirname(file)
    if not os.path.exists(odir): os.makedirs(odir)

    geom = None
    if geom is None:
        try:
            geom = ogr.CreateGeometryFromWkt(wkt)
        except:
            pass

    if geom is None:
        try:
            geom = ogr.CreateGeometryFromJson(wkt)
        except:
            pass

    if geom is not None:
        with open(file, 'w', encoding = 'utf-8') as f: f.write(geom.ExportToJson())
    geom = None

    if os.path.exists(file):
        return(file)
    else:
        return(None)
