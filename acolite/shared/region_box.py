## def region_box
## makes json box of x km around given point
## written by Quinten Vanhellemont, RBINS
## 2021-03-02
## modifications:  2022-04-13 (QV) added return limit option
##

def region_box(site, site_lon, site_lat, box_size = 3, return_limit = False,
                add_box_name = True, sub_dir = None, override=True):
    import acolite as ac
    import os, json

    ## get approximate distance per degree lon/lat
    dlon, dlat = ac.shared.distance_in_ll(site_lat)

    ## make box around station
    lat_off = (box_size/dlat)/2
    lon_off = (box_size/dlon)/2
    limit = [site_lat-lat_off, site_lon-lon_off, site_lat+lat_off, site_lon+lon_off]

    if return_limit:
        return(limit)
    else:
        region_dir = ac.config['data_dir']+'/Regions'
        if sub_dir is not None: region_dir =  '{}/{}'.format(region_dir,sub_dir)
        if not os.path.exists(region_dir): os.makedirs(region_dir)

        jsonregion = '{}-{}x{}km'.format(site,box_size, box_size) if add_box_name else '{}'.format(site)

        jsonf = '{}/{}.geojson'.format(region_dir, jsonregion)
        if (not os.path.exists(jsonf)) or override:
            ## convert limit to linear ring
            coordinates = [[[limit[1],limit[0]],
                            [limit[3],limit[0]],
                            [limit[3],limit[2]],
                            [limit[1],limit[2]],
                            [limit[1],limit[0]] ]]
            data = {"type":"FeatureCollection",
                    "features":[{"type":"Feature","properties":{},
                                 "geometry":{"type":"Polygon","coordinates":coordinates}}]}
            with open(jsonf, 'w') as f: json.dump(data, f)
        return(jsonf)
