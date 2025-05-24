## def bundle_test
## test if given file is EarthCare MSI
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-13
## modifications:

def bundle_test(bundle):
    from netCDF4 import Dataset
    igatts = {}
    with Dataset(bundle) as nc:
        ## get VariableProductHeader
        variables = list(nc.groups['HeaderData']['VariableProductHeader']['MainProductHeader'].variables.keys())
        for k in variables:
            #print(k, nc.groups['HeaderData']['VariableProductHeader']['MainProductHeader'][k][:])
            if k in igatts: print(k)
            igatts[k] = nc.groups['HeaderData']['VariableProductHeader']['MainProductHeader'][k][:]

        ## get SpecificProductHeader
        #variables = list(nc.groups['HeaderData']['VariableProductHeader']['SpecificProductHeader'].variables.keys())
        #for k in variables:
        #    print(k, nc.groups['HeaderData']['VariableProductHeader']['SpecificProductHeader'][k][:])
        #    if k in igatts: print(k)
        #    igatts[k] = nc.groups['HeaderData']['VariableProductHeader']['SpecificProductHeader'][k][:]

        # variables = list(nc.groups['HeaderData']['FixedProductHeader']['Source'].variables.keys())
        # for k in variables:
        #    print(k, nc.groups['HeaderData']['FixedProductHeader']['Source'][k][:])
        #    if k in igatts: print(k)
        #    igatts[k] = nc.groups['HeaderData']['FixedProductHeader']['Source'][k][:]

    ## get mission info
    sensor = None
    if 'description' in igatts:
        if igatts['description'] == 'MSI nominal 1C product':
            platform = 'EarthCare'
            instrument = 'MSI'
            sensor = '{}_{}'.format(platform, instrument)


    return(sensor, igatts)
