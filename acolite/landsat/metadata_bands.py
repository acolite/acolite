## def metadata_bands
## gets landsat band specific metadata
## written by Quinten Vanhellemont, RBINS
## 2021-02-05
## modifications:

def metadata_bands(bundle, meta):
    import os
    fmeta = {}
    if 'PRODUCT_CONTENTS' in meta: ## COLL2
        pk = 'PRODUCT_CONTENTS'
    elif 'PRODUCT_METADATA' in meta: ## COLL1
        pk = 'PRODUCT_METADATA'
    else:
        return(1)

    for k in meta[pk]:
        if 'FILE_NAME' in k:
            fname = meta[pk][k]
            if 'FILE_NAME_' in k:
                par = k[len('_FILE_NAME'):]
            else:
                par = 'BAND_'+k[4:len(k)-len('_FILE_NAME')]
                k = 'FILE_NAME_'+par
            if '.TIF' not in fname: continue
            file = '{}/{}'.format(bundle, fname)
            if os.path.exists(file):
                if 'SOLAR_AZIMUTH' in k:
                    b = "SAA"
                elif 'SOLAR_ZENITH' in k:
                    b = "SZA"
                elif 'SENSOR_AZIMUTH' in k:
                    b = "VAA"
                elif 'SENSOR_ZENITH' in k:
                    b = "VZA"
                else:
                    b = k.split('_')[-1]
                    if 'BAND_' in k:
                        if '_VCID_' in k:
                            b = k[-8:]
                        else:
                            if par != 'BAND_{}'.format(b): continue
                fmeta[b] = {'FILE':file, 'PAR':par}
                for sk in meta.keys():
                    for ssk in meta[sk].keys():
                        if ssk[-len(par):] == par:
                            kb = ssk[0:-len(par)-1]
                            v = meta[sk][ssk]
                            try: v=float(v)
                            except: pass
                            fmeta[b][kb] = v
    return(fmeta)
