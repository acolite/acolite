## def acolite_l1r
## function to convert image from original format to ACOLITE L1R NetCDF
## written by Quinten Vanhellemont, RBINS
## 2021-03-10
## modifications: 2021-12-08 (QV) added nc_projection
##                2021-12-31 (QV) new handling of settings
##                2024-03-14 (QV) update settings handling
##                2025-01-30 (QV) moved polygon limit and limit buffer extension
##                2025-02-04 (QV) fix user and run settings, removed settings from acolite_l1r and l1_convert calls
##                2025-11-04 (QV) added .SAFE and .SEN3 to S2 and S3 input_types

def acolite_l1r(bundle, input_type=None):
    import acolite as ac
    import os, sys, shutil

    ## set up l1r_files list
    l1r_files = []

    ## make bundle a list even if only one is provided
    if type(bundle) != list: bundle = [bundle]
    ac.settings['run']['inputfile'] = bundle
    if 'output' not in ac.settings['run']:
        ac.settings['run']['output'] = os.path.dirname(ac.settings['run']['inputfile'][0])

    ## test path lengths on windows
    if 'win' in sys.platform:
        input_lengths = [len(b) for b in bundle]
        if any([i >= 256 - 82 - 1 for i in input_lengths]): ## 82 chars is the granule band data relative to .SAFE
            print('Warning: Rather long input filename ({} characters)'.format(max(input_lengths)))
            print('This may give issues in Windows due to path length limitations, file(s):')
            for b in bundle: print(len(b), b)

    ## identify bundle
    orig_bundle = [b for b in bundle]
    identification = [ac.acolite.identify_bundle(b, output=ac.settings['run']['output']) for b in bundle]
    input_types = [i[0] for i in identification]
    bundle = [i[1] for i in identification]
    zipped = [i[2] for i in identification]
    extracted_path = [i[3] for i in identification]

    input_type = input_types[0]
    if not all([i == input_type for i in input_types]):
        print('Warning: Multiple input types given: {}'.format(input_types))
    if input_type == None:
        print('{} not recognized.'.format(bundle[0]))

    ## empty setu
    setu = {}

    ################
    ## ACOLITE
    if input_type == 'ACOLITE':
        ## return bundle if ACOLITE type
        l1r_files = bundle
        gatts = ac.shared.nc_gatts(bundle[0])
        setu = ac.acolite.settings.parse(gatts['sensor'], settings = ac.settings['run'])
    ## end ACOLITE
    ################

    ################
    ## Landsat
    if input_type == 'Landsat':
        l1r_files, setu = ac.landsat.l1_convert(bundle)
    ## end Landsat
    ################

    ################
    ## Sentinel-2
    if input_type == 'Sentinel-2 .SAFE':
        l1r_files, setu = ac.sentinel2.l1_convert(bundle)
    ## end Sentinel-2
    ################

    ################
    ## Sentinel-3
    if input_type == 'Sentinel-3 .SEN3':
        l1r_files, setu = ac.sentinel3.l1_convert(bundle)
    if input_type == 'Sentinel-3 .ZARR':
        l1r_files, setu = ac.sentinel3.zarr.l1_convert(bundle)
    ## end Sentinel-3
    ################

    ################
    ## VIIRS
    if input_type == 'VIIRS':
        l1r_files, setu = ac.viirs.l1_convert(bundle)
    ## end VIIRS
    ################

    ################
    ## Pléiades/SPOT
    if input_type == 'Pléiades':
        l1r_files, setu = ac.pleiades.l1_convert(bundle)
    ## end Pléiades/SPOT
    ################

    ################
    ## VENUS
    if input_type == 'VENUS':
        l1r_files, setu = ac.venus.l1_convert(bundle)
    ## end VENUS
    ################

    ################
    ## WorldView
    if input_type == 'WorldView':
        l1r_files, setu = ac.worldview.l1_convert(bundle)
    ## end WorldView
    ################

    ################
    ## CHRIS
    if input_type == 'CHRIS':
        l1r_files, setu = ac.chris.l1_convert(bundle)
    ## end CHRIS
    ################

    ################
    ## PRISMA
    if input_type == 'PRISMA':
        l1r_files, setu = ac.prisma.l1_convert(bundle)
    ## end PRISMA
    ################

    ################
    ## HICO
    if input_type == 'HICO':
        l1r_files, setu = ac.hico.l1_convert(bundle)
    ## end HICO
    ################

    ################
    ## HYPERION
    if input_type == 'HYPERION':
        l1r_files, setu = ac.hyperion.l1_convert(bundle)
    ## end HYPERION
    ################

    ################
    ## DESIS
    if input_type == 'DESIS':
        l1r_files, setu = ac.desis.l1_convert(bundle)
    ## end DESIS
    ################

    ################
    ## PACE
    if input_type == 'PACE':
        l1r_files, setu = ac.pace.l1_convert(bundle)
    ## end PACE
    ################

    ################
    ## SeaDAS L1B
    if input_type == 'SeaDAS':
        l1r_files, setu = ac.seadas.l1_convert(bundle)
    ## end SeaDAS L1B
    ################

    ################
    ## AVHRR
    if input_type == 'AVHRR':
        l1r_files, setu = ac.avhrr.l1_convert(bundle)
    ## end AVHRR
    ################

    ################
    ## Planet
    if input_type == 'Planet':
        l1r_files, setu = ac.planet.l1_convert(bundle)
    ## end Planet
    ################

    ################
    ## GF
    if input_type == 'GF':
        l1r_files, setu = ac.gf.l1_convert(bundle)
    ## end GF
    ################

    ################
    ## AMAZONIA
    if input_type == 'AMAZONIA':
        l1r_files, setu = ac.amazonia.l1_convert(bundle)
    ## end AMAZONIA
    ################

    ################
    ## FORMOSAT
    if input_type == 'FORMOSAT':
        l1r_files, setu = ac.formosat.l1_convert(bundle)
    ## end FORMOSAT
    ################

    ################
    ## SDGSAT
    if input_type == 'SDGSAT':
        l1r_files, setu = ac.sdgsat.l1_convert(bundle)
    ## end SDGSAT
    ################

    ################
    ## IKONOS
    if input_type == 'IKONOS':
        l1r_files, setu = ac.ikonos.l1_convert(bundle)
    ## end IKONOS
    ################

    ################
    ## DEIMOS
    if input_type == 'DEIMOS':
        l1r_files, setu = ac.deimos.l1_convert(bundle)
    ## end DEIMOS
    ################

    ################
    ## ECOSTRESS
    if input_type == 'ECOSTRESS':
        l1r_files, setu = ac.ecostress.l1_convert(bundle)
    ## end ECOSTRESS
    ################

    ################
    ## ENMAP
    if input_type == 'ENMAP':
        l1r_files, setu = ac.enmap.l1_convert(bundle)
    ## end ENMAP
    ################

    ################
    ## EMIT
    if input_type == 'EMIT':
        l1r_files, setu = ac.emit.l1_convert(bundle)
    ## end EMIT
    ################

    ################
    ## Tanager
    if input_type == 'Tanager':
        l1r_files, setu = ac.tanager.l1_convert(bundle)
    ## end Tanager
    ################

    ################
    ## EarthCare
    if input_type == 'EarthCare':
        l1r_files, setu = ac.earthcare.l1_convert(bundle)
    ## end EarthCare
    ################

    ################
    ## DIMAP
    if input_type == 'DIMAP':
        l1r_files, setu = ac.dimap.l1_convert(bundle)
    ## end DIMAP
    ################

    ################
    ## S2Resampling
    if input_type == 'S2Resampling':
        l1r_files, setu = ac.s2resampling.l1_convert(bundle)
    ## end S2Resampling
    ################

    ################
    ## SEVIRI
    if input_type == 'SEVIRI':
        l1r_files, setu = ac.seviri.l1_convert(bundle)
    ## end SEVIRI
    ################

    ################
    ## FCI
    if input_type == 'FCI':
        l1r_files, setu = ac.fci.l1_convert(bundle)
    ## end FCI
    ################

    ################
    ## AHI
    if input_type == 'Himawari':
        l1r_files, setu = ac.himawari.l1_convert(bundle)
    ## end AHI
    ################

    ################
    ## ABI
    if input_type == 'GOES':
        l1r_files, setu = ac.goes.l1_convert(bundle)
    ## end ABI
    ################

    ################
    ## GOCI
    if input_type == 'GOCI':
        l1r_files, setu = ac.goci.l1_convert(bundle)
    ## end GOCI
    ################

    ################
    ## Haiyang
    if input_type == 'HAIYANG':
        l1r_files, setu = ac.haiyang.l1_convert(bundle)
    ## end Haiyang
    ################

    ################
    ## Huanjing
    if input_type == 'HUANJING':
        l1r_files, setu = ac.huanjing.l1_convert(bundle)
    ## end Huanjing
    ################

    ################
    ## HYPSO
    if input_type == 'HYPSO':
        l1r_files, setu = ac.hypso.l1_convert(bundle)
    ## end HYPSO
    ################

    ################
    ## Wyvern
    if input_type == 'Wyvern':
        l1r_files, setu = ac.wyvern.l1_convert(bundle)
    ## end Wyvern
    ################

    ################
    ## OpenCosmos
    if input_type == 'OpenCosmos':
        l1r_files, setu = ac.opencosmos.l1_convert(bundle)
    ## end OpenCosmos
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
