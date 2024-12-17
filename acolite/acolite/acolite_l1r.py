## def acolite_l1r
## function to convert image from original format to ACOLITE L1R NetCDF
## written by Quinten Vanhellemont, RBINS
## 2021-03-10
## modifications: 2021-12-08 (QV) added nc_projection
##                2021-12-31 (QV) new handling of settings
##                2024-03-14 (QV) update settings handling

def acolite_l1r(bundle, settings = None, input_type=None):
    import acolite as ac
    import os, sys, shutil

    ## parse run settings
    if settings is not None:
        ac.settings['user'] = ac.acolite.settings.parse(None, settings=settings, merge=False)
        for k in ac.settings['user']: ac.settings['run'][k] = ac.settings['user'][k]
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## set up l1r_files list
    l1r_files = []

    ## make bundle a list even if only one is provided
    if type(bundle) != list: bundle = [bundle]
    setu['inputfile'] = bundle

    ## test path lengths on windows
    if 'win' in sys.platform:
        input_lengths = [len(b) for b in bundle]
        if any([i >= 256 - 82 - 1 for i in input_lengths]): ## 82 chars is the granule band data relative to .SAFE
            print('Warning: Rather long input filename ({} characters)'.format(max(input_lengths)))
            print('This may give issues in Windows due to path length limitations, file(s):')
            for b in bundle: print(len(b), b)

    ## set output directory
    if 'output' not in setu:
        setu['output'] = os.path.dirname(setu['inputfile'][0])

    ## identify bundle
    orig_bundle = [b for b in bundle]
    identification = [ac.acolite.identify_bundle(b, output=setu['output']) for b in bundle]
    input_types = [i[0] for i in identification]
    bundle = [i[1] for i in identification]
    zipped = [i[2] for i in identification]
    extracted_path = [i[3] for i in identification]

    input_type = input_types[0]
    if not all([i == input_type for i in input_types]):
        print('Warning: Multiple input types given: {}'.format(input_types))
    if input_type == None:
        print('{} not recognized.'.format(bundle[0]))

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
    ## VIIRS
    if input_type == 'VIIRS':
        l1r_files, setu = ac.viirs.l1_convert(bundle, settings = setu)
    ## end VIIRS
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
    ## PACE
    if input_type == 'PACE':
        l1r_files, setu = ac.pace.l1_convert(bundle, settings = setu)
    ## end PACE
    ################

    ################
    ## SeaDAS L1B
    if input_type == 'SeaDAS':
        l1r_files, setu = ac.seadas.l1_convert(bundle, settings = setu)
    ## end SeaDAS L1B
    ################

    ################
    ## AVHRR
    if input_type == 'AVHRR':
        l1r_files, setu = ac.avhrr.l1_convert(bundle, settings = setu)
    ## end AVHRR
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

    ################
    ## SDGSAT
    if input_type == 'SDGSAT':
        l1r_files, setu = ac.sdgsat.l1_convert(bundle, settings = setu)
    ## end SDGSAT
    ################

    ################
    ## IKONOS
    if input_type == 'IKONOS':
        l1r_files, setu = ac.ikonos.l1_convert(bundle, settings = setu)
    ## end IKONOS
    ################

    ################
    ## ECOSTRESS
    if input_type == 'ECOSTRESS':
        l1r_files, setu = ac.ecostress.l1_convert(bundle, settings = setu)
    ## end ECOSTRESS
    ################

    ################
    ## ENMAP
    if input_type == 'ENMAP':
        l1r_files, setu = ac.enmap.l1_convert(bundle, settings = setu)
    ## end ENMAP
    ################

    ################
    ## EMIT
    if input_type == 'EMIT':
        l1r_files, setu = ac.emit.l1_convert(bundle, settings = setu)
    ## end EMIT
    ################

    ################
    ## DIMAP
    if input_type == 'DIMAP':
        l1r_files, setu = ac.dimap.l1_convert(bundle, settings = setu)
    ## end DIMAP
    ################

    ################
    ## S2Resampling
    if input_type == 'S2Resampling':
        l1r_files, setu = ac.s2resampling.l1_convert(bundle, settings = setu)
    ## end S2Resampling
    ################

    ################
    ## SEVIRI
    if input_type == 'SEVIRI':
        l1r_files, setu = ac.seviri.l1_convert(bundle, settings = setu)
    ## end SEVIRI
    ################

    ################
    ## Haiyang
    if input_type == 'HAIYANG':
        l1r_files, setu = ac.haiyang.l1_convert(bundle, settings = setu)
    ## end Haiyang
    ################

    ################
    ## HYPSO
    if input_type == 'HYPSO':
        l1r_files, setu = ac.hypso.l1_convert(bundle, settings = setu)
    ## end HYPSO
    ################

    ## remove extracted files
    for i, im in enumerate(identification):
        try:
            if (im[2]) & (setu['delete_extracted_input']):
                ob = os.path.abspath(orig_bundle[i])
                eb = os.path.abspath(extracted_path[i])
                if ob != eb:
                    if os.path.exists(eb) & (extracted_path[i] not in orig_bundle):
                        print('Deleting {}'.format(eb))
                        shutil.rmtree(eb)
        except KeyError:
            print('Not deleting extracted file as "delete_extracted_input" is not in settings.')
        except BaseException as err:
            print("Error removing extracted files {}, {}".format(err, type(err)))

    return(l1r_files, setu, bundle)
