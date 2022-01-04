## def acolite_l1r
## function to convert image from original format to ACOLITE L1R NetCDF
## written by Quinten Vanhellemont, RBINS
## 2021-03-10
## modifications: 2021-12-08 (QV) added nc_projection
##                2021-12-31 (QV) new handling of settings

def acolite_l1r(bundle, setu, input_type=None):
    import acolite as ac
    import os

    ## set up l1r_files list
    l1r_files = []

    ## identify bundle
    if type(bundle) != list: bundle = [bundle]
    bundle.sort()
    input_type = ac.acolite.identify_bundle(bundle[0])
    if input_type is None:
        print('{} not recognized.'.format(bundle[0]))

    setu['inputfile'] = bundle
    ## set output directory
    if 'output' not in setu:
        setu['output'] = os.path.dirname(setu['inputfile'][0])

    ################
    ## ACOLITE
    if input_type == 'ACOLITE':
        ## return bundle if ACOLITE type
        l1r_files = bundle
        gatts = ac.shared.nc_gatts(bundle[0])
        setu = ac.acolite.settings.parse(gatts['sensor'], settings = setu)
    ## end ACOLITE
    ################

    ################
    ## Landsat
    if input_type == 'Landsat':
        l1r_files, setu = ac.landsat.l1_convert(bundle, settings = setu)
    ## end Landsat
    ################

    ################
    ## Sentinel-2
    if input_type == 'Sentinel-2':
        l1r_files, setu = ac.sentinel2.l1_convert(bundle, settings = setu)
    ## end Sentinel-2
    ################

    ################
    ## Sentinel-3
    if input_type == 'Sentinel-3':
        l1r_files, setu = ac.sentinel3.l1_convert(bundle, settings = setu)
    ## end Sentinel-3
    ################

    ################
    ## Pléiades/SPOT
    if input_type == 'Pléiades':
        l1r_files, setu = ac.pleiades.l1_convert(bundle, settings = setu)
    ## end Pléiades/SPOT
    ################

    ################
    ## VENUS
    if input_type == 'VENUS':
        l1r_files, setu = ac.venus.l1_convert(bundle, settings = setu)
    ## end VENUS
    ################

    ################
    ## WorldView
    if input_type == 'WorldView':
        if 'inputfile_swir' not in setu: setu['inputfile_swir'] = None
        l1r_files, setu = ac.worldview.l1_convert(bundle, inputfile_swir = setu['inputfile_swir'], settings = setu)
    ## end WorldView
    ################

    ################
    ## CHRIS
    if input_type == 'CHRIS':
        l1r_files, setu = ac.chris.l1_convert(bundle, settings = setu)
    ## end CHRIS
    ################

    ################
    ## PRISMA
    if input_type == 'PRISMA':
        l1r_files, setu = ac.prisma.l1_convert(bundle, settings = setu)
    ## end PRISMA
    ################

    ################
    ## HICO
    if input_type == 'HICO':
        l1r_files, setu = ac.hico.l1_convert(bundle, settings = setu)
    ## end HICO
    ################

    ################
    ## HYPERION
    if input_type == 'HYPERION':
        l1r_files, setu = ac.hyperion.l1_convert(bundle, settings = setu)
    ## end HYPERION
    ################

    ################
    ## DESIS
    if input_type == 'DESIS':
        l1r_files, setu = ac.desis.l1_convert(bundle, settings = setu)
    ## end DESIS
    ################

    ################
    ## Planet
    if input_type == 'Planet':
        l1r_files, setu = ac.planet.l1_convert(bundle, settings = setu)
    ## end Planet
    ################

    ################
    ## GF
    if input_type == 'GF':
        l1r_files, setu = ac.gf.l1_convert(bundle, settings = setu)
    ## end GF
    ################


    return(l1r_files, setu)
