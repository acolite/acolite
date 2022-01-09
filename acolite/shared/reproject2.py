## def reproject2
## reproject data from one lon/lat combination to another
## written by Quinten Vanhellemont, RBINS
## QV 2017-07-17
## modifications: 2022-01-09 (QV) added radius_of_influence as keyword

def reproject2(data, lon0, lat0, lon1, lat1, fill=-9999, nearest=True, radius_of_influence=100):
    from pyresample import image, geometry

    source_def = geometry.SwathDefinition(lons=lon0,lats=lat0)
    target_def = geometry.SwathDefinition(lons=lon1,lats=lat1)

    if nearest:
        source_con = image.ImageContainerNearest(data,source_def,
                                                radius_of_influence=radius_of_influence)
        target_con = source_con.resample(target_def)
        target_con.fill_value=fill
        result = target_con.image_data
    else:
        from pyresample import kd_tree
        result = kd_tree.resample_gauss(source_def, data, target_def,
                                        radius_of_influence=radius_of_influence,
                                        neighbours=1,sigmas=0, fill_value=fill)

    return result
