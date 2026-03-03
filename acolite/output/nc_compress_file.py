## def nc_compress_file
## applies compression to NetCDF file that was previously written
## written by Quinten Vanhellemont, RBINS
## 2026-03-03
## modifications:

def nc_compress_file(inputfile, output = None, replace = False, delete = True,
                    temp_ext = '_temp.nc', comp_ext = '_comp.nc',
                    netcdf_compression_level = None, netcdf_compression_least_significant_digit = None):
    import os
    import acolite as ac

    ## override compression setting
    ac.settings['run']['netcdf_compression'] = True
    if netcdf_compression_level is not None:
        ac.settings['run']['netcdf_compression_level'] = netcdf_compression_level
    if netcdf_compression_least_significant_digit is not None:
        ac.settings['run']['netcdf_compression_least_significant_digit'] = netcdf_compression_least_significant_digit

    ## make new outputile
    bn = os.path.basename(inputfile)
    bn_, ext_ = os.path.splitext(bn)

    ## target directory - used if replace is not set
    if output is None:
        dn = os.path.dirname(inputfile)
    else:
        dn = output

    ## do we make a new file or replace the existing file
    if replace:
        dn_ = os.path.dirname(inputfile)
        tmpfile = '{}/{}{}'.format(dn_, bn_, temp_ext)
        os.rename(inputfile, tmpfile)
        ofile = '{}'.format(inputfile)
        inputfile = '{}'.format(tmpfile)
    else:
        tmpfile = None
        ofile = '{}/{}{}'.format(dn, bn_, comp_ext)

    ## open inputfile
    gem = ac.gem.gem(inputfile)

    ## create dataset
    gemo = ac.gem.gem(ofile, new = True)

    ## first copy all input data
    gemo.gatts = {k: gem.gatts[k] for k in gem.gatts} ## set output attributes
    gemo.gatts['ofile'] = ofile
    if gem.nc_projection is not None: gemo.nc_projection = gem.nc_projection

    for ds in gem.datasets:
        if gem.nc_projection is not None:
            if ds in gem.nc_projection_keys: continue

        try:
            d, da = gem.data(ds, attributes = True)
        except:
            print('Could not read dataset {}'.format(ds))
            continue

        ## output dataset
        gemo.write(ds, d, ds_att = da)
        if ac.settings['run']['verbosity'] > 2:
            print('Wrote {} ({}) to {}'.format(ds, d.shape, ofile))
    gem.close()

    if (replace) & (delete):
        if os.path.exists(tmpfile):
            os.remove(tmpfile)

    return(ofile)
