## def acolite_l1r
## function to convert image from original format to ACOLITE L1R NetCDF
## written by Quinten Vanhellemont, RBINS
## 2021-03-10
## modifications: 2021-12-08 (QV) added nc_projection
##                2021-12-31 (QV) new handling of settings

def acolite_l1r(bundle, setu, input_type=None):
    import acolite as ac
    import os, sys

    ## set up l1r_files list
    l1r_files = []

    ## make bundle a list even if only one is provided
    if type(bundle) != list: bundle = [bundle]
    bundle.sort()
    setu['inputfile'] = bundle

    ## test path lengths on windows
    if 'win' in sys.platform:
        input_lengths = [len(b) for b in bundle]
        if any([i >= 256 - 82 for i in input_lengths]): ## 82 chars is the granule band data relative to .SAFE
            print('Warning: Rather long input filename ({} characters)'.format(max(input_lengths)))
            print('This may give issues in Windows due to path length limitations, file(s):')
            for b in bundle: print(len(b), b)

    ## identify bundle
    input_types = [ac.acolite.identify_bundle(b) for b in bundle]
    input_type = input_types[0]
    if not all([i == input_type for i in input_types]):
        print('Warning: Multiple input types given: {}'.format(input_types))
    if input_type == None:
        print('{} not recognized.'.format(bundle[0]))

    ## set output directory
    if 'output' not in setu:
        setu['output'] = os.path.dirname(setu['inputfile'][0])

    if 'limit' in setu:
        if setu['limit'] is not None:
            if len(setu['limit']) != 4:
                print('ROI limit should be four elements in decimal degrees: limit=S,W,N,E')
                print('Provided in the settings:', setu['limit'])
                return()

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

    ################
    ## AMAZONIA
    if input_type == 'AMAZONIA':
        l1r_files, setu = ac.amazonia.l1_convert(bundle, settings = setu)
    ## end AMAZONIA
    ################

    ################
    ## FORMOSAT
    if input_type == 'FORMOSAT':
        l1r_files, setu = ac.formosat.l1_convert(bundle, settings = setu)
    ## end FORMOSAT
    ################

    return(l1r_files, setu)
