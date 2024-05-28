## def auxiliary_json
## reads auxiliary data included from S2 PB004 onward and stores as json file
## written by Quinten Vanhellemont, RBINS
## 2024-05-28
## modifications:

def auxiliary_json(bundle, output = None):
    import os, json
    import acolite as ac

    ## parse bundle and get granule
    safe_files = ac.sentinel2.safe_test(bundle)
    if 'granules' not in safe_files:
        print('{} not recognised'.format(bundle))
        return

    if output is None:
        odir = os.path.dirname(granule)
    else:
        odir = '{}'.format(output)

    ## only one granule in L1C
    granule = safe_files['granules'][0]

    ## get granule auxillary data
    data = ac.sentinel2.auxiliary(bundle, granule)

    ## extract grid points
    winds = ((data['u10']['values']**2 + data['v10']['values']**2) ** 0.5).flatten()
    aots = (data['aod550']['values']).flatten()
    lats =  (data['aod550']['latitudes']).flatten()
    lons =  (data['aod550']['longitudes']).flatten()

    ## write parameter JSON file
    parameters = {'wind': [float(v) for v in winds],
                  'lon': [float(v) for v in lons],
                  'lat': [float(v) for v in lats],
                  'aot': [float(v) for v in aots]}

    ## output json
    ofile_json = '{}/{}_AUX_TPG.json'.format(odir, os.path.basename(os.path.splitext(bundle)[0]))
    if not os.path.exists(odir): os.makedirs(odir)
    with open(ofile_json, 'w') as f:
        json.dump(parameters, f)

    print('Wrote {}'.format(ofile_json))
    return(ofile_json)
