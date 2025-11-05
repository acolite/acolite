## def meta_parse
## parse metadata from from zarr meta
## written by Quinten Vanhellemont, RBINS
## 2025-11-04
## modifications:

def meta_parse(z):
    import acolite as ac

    ## set up metadata dict
    meta = {'sensor': None}

    ## parse zarr metadata
    if type(z) is not dict:
        md = ac.zarr.meta(z)
    else:
        md = {k: z[k] for k in z}

    ## extract platform information
    prop = md['attributes']['stac_discovery']['properties']
    if prop['platform'].startswith('sentinel'):
        meta['platform'] = 'S'+prop['platform'][-2:].upper()
        if 'instrument' in prop:
            meta['instrument'] = prop['instrument'].upper()
        elif 'instruments' in prop:
            meta['instrument'] = prop['instruments'][0].upper()
        else:
            print('Instrument could not be determined.')
            print(prop.keys())
        meta['sensor'] = '{}_{}'.format(meta['platform'], meta['instrument'])

    ## identifier
    for k in ['id', 'collection']:
        if k in md['attributes']['stac_discovery']:
            meta[k] = md['attributes']['stac_discovery'][k]

    ## extract keys
    for k in ['datetime', 'start_datetime', 'end_datetime',
              'sat:absolute_orbit', 'sat:relative_orbit', 'sat:orbit_state',
              'proj:shape', 'processing:level', 'product:type']:
        if k in prop: meta[k] = prop[k]

    ## extract other metadata
    other = md['attributes']['other_metadata']
    for k in ['data_information', 'absolute_pass_number', 'relative_pass_number', 'cycle_number'
              'UTM_zone_identification', 'horizontal_CRS_code', 'horizontal_CRS_name',
              'mean_sun_azimuth_angle_in_deg_for_all_bands_all_detectors', 'mean_sun_zenith_angle_in_deg_for_all_bands_all_detectors',
              'reflectance_correction_factor_from_the_Sun_Earth_distance_variation_computed_using_the_acquisition_date']:
        if k in other: meta[k] = other[k]

    ## convert lists of len 1
    for k in meta:
        if type(meta[k]) == list:
            if len(meta[k]) == 1: meta[k] = meta[k][0]

    return(meta)
