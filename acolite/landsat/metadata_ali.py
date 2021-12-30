## def metadata_ali
## updates metadata for EO1/ALI
## written by Quinten Vanhellemont, RBINS
## 2021-12-24
## modifications:

def metadata_ali(meta):
    import acolite as ac
    import datetime, dateutil.parser

    if 'PRODUCT_METADATA' in meta: ## COLL1
        pk = 'PRODUCT_METADATA'
        ik = 'IMAGE_ATTRIBUTES'
        rk = 'PRODUCT_METADATA'

        if (meta[pk]['SPACECRAFT_ID'] == 'EO1') & (meta[pk]['SENSOR_ID'] == 'ALI'):
            meta['PRODUCT_METADATA']['DATA_TYPE'] = meta['PRODUCT_METADATA']['PRODUCT_TYPE']
            meta[pk]['WRS_PATH'] = ''
            meta[pk]['WRS_ROW'] = ''
            meta[pk]['DATE_ACQUIRED'] = meta[pk]['ACQUISITION_DATE']
            st = dateutil.parser.parse(meta[pk]['DATE_ACQUIRED']+'T'+meta[pk]['START_TIME'].split(' ')[-1])
            et = dateutil.parser.parse(meta[pk]['DATE_ACQUIRED']+'T'+meta[pk]['END_TIME'].split(' ')[-1])
            stime = st + datetime.timedelta(seconds=(et-st).seconds/2)
            meta[pk]['SCENE_CENTER_TIME'] = stime.strftime('%H:%M:%S')
            meta[rk]['REFLECTIVE_SAMPLES'] = meta[rk]['PRODUCT_SAMPLES_REF']
            meta[rk]['REFLECTIVE_LINES'] = meta[rk]['PRODUCT_LINES_REF']

            meta['PROJECTION_PARAMETERS']['GRID_CELL_SIZE_REFLECTIVE'] = meta['PROJECTION_PARAMETERS']['GRID_CELL_SIZE_REF']
            meta['PROJECTION_PARAMETERS']['UTM_ZONE'] = meta['UTM_PARAMETERS']['ZONE_NUMBER']

            spos = ac.shared.sun_position(stime, 0, 0)
            meta[ik] = {'SUN_ELEVATION':  meta['PRODUCT_PARAMETERS']['SUN_ELEVATION'],
                                'SUN_AZIMUTH':  meta['PRODUCT_PARAMETERS']['SUN_AZIMUTH'],
                                'EARTH_SUN_DISTANCE': spos['distance']}

            ## add band info
            f0_b = {'1': 17237.94933749016,
                    '2': 18562.143292474477,
                    '3': 19947.039214591663,
                    '4': 18068.91057456214,
                    '5': 15362.385878114905,
                    '6': 11443.287934895432,
                    '7': 9548.229474857062,
                    '8': 4524.315423791478,
                    '9': 2350.783374529712,
                    '10': 823.8106067793342}

            meta['ALI_BANDS'] = {}
            for b in range(1,11):
                meta['ALI_BANDS']['RADIANCE_MULT_BAND_{}'.format(b)] = meta['RADIANCE_SCALING']['BAND{}_SCALING_FACTOR'.format(b)]
                meta['ALI_BANDS']['RADIANCE_ADD_BAND_{}'.format(b)] = meta['RADIANCE_SCALING']['BAND{}_OFFSET'.format(b)]

                meta['ALI_BANDS']['QUANTIZE_CAL_MIN_BAND_{}'.format(b)] = 1
                meta['ALI_BANDS']['QUANTIZE_CAL_MAX_BAND_{}'.format(b)] = 32767

                meta['ALI_BANDS']['se_distance_BAND_{}'.format(b)] = spos['distance']
                meta['ALI_BANDS']['f0_BAND_{}'.format(b)] = f0_b['{}'.format(b)]/10

    return(meta)
