## def acolite_l1r
## function to convert image from original format to ACOLITE L1R NetCDF
## written by Quinten Vanhellemont, RBINS
## 2021-03-10
## modifications:
##

def acolite_l1r(bundle, output, setu, input_type=None):
    import acolite as ac
    import os

    ## set up l1r_files list
    l1r_files = []

    ## identify bundle
    if type(bundle) != list: bundle = [bundle]
    bundle.sort()
    input_type = ac.acolite.identify_bundle(bundle[0])
    if input_type is None:
        if setu['verbosity'] > 0: print('{} not recognized.'.format(bundle[0]))

    ## set output directory
    if output is None:
        output_ = setu['output'] if setu['output'] is not None else os.path.dirname(bundle[0])
    else:
        output_ = output

    print(output_)
    ################
    ## ACOLITE
    if input_type == 'ACOLITE':
        ## return bundle if ACOLITE type
        l1r_files = bundle
    ## end ACOLITE
    ################

    ################
    ## Landsat
    if input_type == 'Landsat':
        l1r_files = ac.landsat.l1_convert(bundle, output=output_,
                                    limit=setu['limit'], poly=setu['polygon'],
                                    merge_tiles = setu['merge_tiles'],
                                    merge_zones = setu['merge_zones'],
                                    extend_region = setu['extend_region'],
                                    output_geolocation=setu['output_geolocation'],
                                    output_xy=setu['output_xy'],
                                    output_geometry=setu['output_geometry'],
                                    vname = setu['region_name'],
                                    verbosity=setu['verbosity'])
    ## end Landsat
    ################

    ################
    ## Sentinel-2
    if input_type == 'Sentinel-2':
        l1r_files = ac.sentinel2.l1_convert(bundle, output=output_,
                                     s2_target_res=setu['s2_target_res'],
                                     geometry_type=setu['geometry_type'],
                                     geometry_res=setu['geometry_res'],
                                     limit=setu['limit'], poly=setu['polygon'],
                                     merge_tiles = setu['merge_tiles'],
                                     merge_zones = setu['merge_zones'],
                                     extend_region = setu['extend_region'],
                                     output_geolocation=setu['output_geolocation'],
                                     output_xy=setu['output_xy'],
                                     output_geometry=setu['output_geometry'],
                                     vname = setu['region_name'],
                                     verbosity=setu['verbosity'])
    ## end Sentinel-2
    ################

    ################
    ## Sentinel-3
    if input_type == 'Sentinel-3':
        l1r_files = ac.sentinel3.l1_convert(bundle, output=output_,
                                            smile_correction = setu['smile_correction'],
                                            use_tpg = setu['use_tpg'],
                                            limit=setu['limit'], poly=setu['polygon'],
                                            output_geolocation=setu['output_geolocation'],
                                            output_geometry=setu['output_geometry'],
                                            vname = setu['region_name'],
                                            verbosity=setu['verbosity'])
    ## end Sentinel-3
    ################

    ################
    ## Pléiades/SPOT
    if input_type == 'Pléiades':
        l1r_files = ac.pleiades.l1_convert(bundle, output=output_,
                                           limit=setu['limit'], poly=setu['polygon'],
                                           output_geolocation=setu['output_geolocation'],
                                           vname = setu['region_name'],
                                           verbosity=setu['verbosity'])
    ## end Pléiades/SPOT
    ################

    ################
    ## WorldView
    if input_type == 'WorldView':
        l1r_files = ac.worldview.l1_convert(bundle, output=output_,
                                            inputfile_swir = setu['inputfile_swir'],
                                            limit=setu['limit'], poly=setu['polygon'],
                                            output_geolocation=setu['output_geolocation'],
                                            vname = setu['region_name'],
                                            verbosity=setu['verbosity'])
    ## end WorldView
    ################

    ################
    ## Planet
    if input_type == 'Planet':
        l1r_files = ac.planet.l1_convert(bundle, output=output_,
                                         limit=setu['limit'], poly=setu['polygon'],
                                         merge_tiles = setu['merge_tiles'],
                                         merge_zones = setu['merge_zones'],
                                         extend_region = setu['extend_region'],
                                         output_geolocation=setu['output_geolocation'],
                                         output_xy=setu['output_xy'],
                                         vname = setu['region_name'],
                                         verbosity=setu['verbosity'])
    ## end Planet
    ################

    return(l1r_files)
