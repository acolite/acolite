## def identify_bundle
## function to identify image that needs to be processed
## written by Quinten Vanhellemont, RBINS
## 2021-03-10
## modifications: 2021-04-08 (QV) added VENUS
##                2022-07-21 (QV) added check for tar and zip files
##                2022-10-26 (QV) update Planet multi-scene zip file handling

def identify_bundle(bundle, input_type = None, output = None):
    import os, glob, shutil, zipfile
    import acolite as ac

    zipped = False
    orig_bundle = '{}'.format(bundle)
    extracted_path = None

    while input_type is None:
        if not os.path.exists(bundle):
            print('Input file {} does not exist'.format(bundle))
            break ## exit loop if path does not exist

        ## test if zip/tar file
        bn, ext = os.path.splitext(bundle)
        if os.path.isfile(bundle) & (ext.lower() in ['.zip', '.tar', '.gz', '.tgz', '.bz2']):
            targ_bundle, extracted_path = ac.shared.extract_bundle(bundle, output=output, verbosity=2)
            if targ_bundle is not None:
                print(targ_bundle)
                bundle = '{}'.format(targ_bundle)
                zipped = True

        ################
        ## ACOLITE
        try:
            gatts = ac.shared.nc_gatts(bundle)
            datasets = ac.shared.nc_datasets(bundle)
            rhot_ds = [ds for ds in datasets if 'rhot_' in ds]
            rhot_ds += [ds for ds in datasets if ds[0:2]=='bt']
            if (gatts['generated_by'] == 'ACOLITE') & (len(rhot_ds) != 0):
                input_type = 'ACOLITE'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end ACOLITE
        ################

        ################
        ## Landsat
        try:
            mtl = glob.glob('{}/{}'.format(bundle, '*MTL.txt'))
            mtl += glob.glob('{}/{}'.format(bundle, '*MTL_L1T.TXT'))
            mtl += glob.glob('{}/{}'.format(bundle, '*MTL_L1GST.TXT'))
            if len(mtl)>0:
                meta = ac.landsat.metadata_read(mtl[0])
                ## get relevant data from meta
                if 'PRODUCT_CONTENTS' in meta: pk = 'IMAGE_ATTRIBUTES'## COLL2
                elif 'PRODUCT_METADATA' in meta: pk = 'PRODUCT_METADATA'## COLL1
                spacecraft_id, sensor_id = meta[pk]['SPACECRAFT_ID'],meta[pk]['SENSOR_ID']
                if (spacecraft_id in ['LANDSAT_5', 'LANDSAT_7', 'LANDSAT_8', 'LANDSAT_9']) | ((spacecraft_id == 'EO1') & (sensor_id == 'ALI')):
                    input_type = 'Landsat'
                    break ## exit loop
                if (spacecraft_id in ['LANDSAT_1', 'LANDSAT_2', 'LANDSAT_3', 'LANDSAT_4', 'LANDSAT_5']) & (sensor_id == 'MSS'):
                    input_type = 'Landsat'
        except:
            pass ## continue to next sensor
        ## end Landsat
        ################

        ################
        ## Sentinel-2
        try:
            safe_files = ac.sentinel2.safe_test(bundle)
            granule = safe_files['granules'][0]
            meta, band_data= ac.sentinel2.metadata_scene(safe_files['metadata']['path'])
            if meta['SPACECRAFT_NAME'] in ['Sentinel-2A', 'Sentinel-2B', 'Sentinel-2C']:
                input_type = 'Sentinel-2'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end Sentinel-2
        ################

        ################
        ## Sentinel-3
        try:
            if zipped & (len(glob.glob('{}/*.SEN3'.format(bundle))) == 1):
                dfiles = glob.glob('{}/*.SEN3/*.nc'.format(bundle))
                dfiles.sort()
            else:
                dfiles = glob.glob('{}/*.nc'.format(bundle))
                dfiles.sort()
            gatts = ac.shared.nc_gatts(dfiles[0])
            if 'OLCI Level 1b Product' in gatts['title']:
                input_type = 'Sentinel-3'
                break ## exit loop
            elif 'MERIS Level 1b Product' in gatts['title']:
                input_type = 'Sentinel-3'
                break ## exit loop
            elif 'S3 SLSTR L1' in gatts['title']:
                input_type = 'SLSTR'
                input_type = 'Sentinel-3'

                break ## exit loop
            else:
                print(gatts['title'])
        except:
            pass ## continue to next sensor
        ## end Sentinel-3
        ################

        ################
        ## VIIRS
        try:
            ret = ac.viirs.bundle_test(bundle)
            if (ret is not None):
                input_type = 'VIIRS'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end VIIRS
        ################

        ################
        ## S2Resampling
        try:
            ret = ac.s2resampling.bundle_test(bundle)
            if (ret is not None):
                input_type = 'S2Resampling'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end S2Resampling
        ################

        ################
        ## Pléiades/SPOT
        try:
            ifiles,mfiles,pifiles,pmfiles = ac.pleiades.bundle_test(bundle, listpan=True)
            mfiles_set = set(mfiles)
            for mfile in mfiles_set: meta = ac.pleiades.metadata_parse(mfile)
            if meta['satellite'] in ['Pléiades', 'SPOT']:
                input_type = 'Pléiades'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end Pléiades/SPOT
        ################

        ################
        ## VENUS
        try:
            meta = ac.venus.metadata_parse(bundle)
            if meta['PLATFORM'] == 'VENUS':
                input_type = 'VENUS'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end VENUS
        ################

        ################
        ## WorldView
        try:
            metafile = ac.worldview.bundle_test(bundle)
            meta = ac.worldview.metadata_parse(metafile)
            if meta['satellite'] in ['WorldView2', 'WorldView3', 'QuickBird2', 'GeoEye1',\
                                     'LG01', 'LG02', 'LG03', 'LG04', 'LG05', 'LG06', ]:
                input_type = 'WorldView'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end WorldView
        ################

        ################
        ## CHRIS
        try:
            gains, mode_info = ac.chris.vdata(bundle)
            input_type = 'CHRIS'
            break ## exit loop
        except:
            pass ## continue to next sensor
        ## end CHRIS
        ################

        ################
        ## PRISMA
        try:
            gatts = ac.prisma.attributes(bundle)
            if gatts['Product_ID'] == b'PRS_L1_STD':
                input_type = 'PRISMA'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end PRISMA
        ################

        ################
        ## HICO
        try:
            gatts = ac.hico.attributes(bundle)
            if gatts['Instrument_Short_Name'] == 'hico':
                input_type = 'HICO'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end HICO
        ################

        ################
        ## HYPERION
        try:
            metadata = ac.hyperion.metadata(bundle)
            if metadata['PRODUCT_METADATA']['SENSOR_ID'] == 'HYPERION':
                input_type = 'HYPERION'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end HYPERION
        ################

        ################
        ## DESIS
        try:
            metafile, imagefile = ac.desis.bundle_test(bundle)
            #headerfile = imagefile.replace('.tif', '.hdr')
            # header = ac.desis.hdr(headerfile)
            meta = ac.desis.metadata(metafile)
            if (meta['mission'] == 'DESIS') & (meta['productType'] in ['L1B', 'L1C']):
                input_type = 'DESIS'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end DESIS
        ################

        ################
        ## PACE
        try:
            gatts = ac.shared.nc_gatts(bundle)
            if 'PACE' in gatts['title']: platform = 'PACE'
            if 'OCI Level-2 Data' in gatts['title']: platform = 'PACE'
            if '{}_{}'.format(platform, gatts['instrument']) == 'PACE_OCI':
                input_type = 'PACE'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end PACE
        ################

        ################
        ## SeaDAS L1B
        try:
            gatts = ac.shared.nc_gatts(bundle)
            if (gatts['processing_level'] == 'L1B') & \
                (gatts['instrument'] in ['HAWKEYE']):
                input_type = 'SeaDAS'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end SeaDAS L1B
        ################

        ################
        ## Planet data
        ## unzip files if needed
        try:
            ## test files
            files = ac.planet.bundle_test(bundle)
            ## files now contains scene_id keys
            flist = list(files.keys())
            flist.sort()
            fk = flist[0]
            if 'metadata' in files[fk]:
                metafile = files[fk]['metadata']['path']
            elif 'metadata_json' in files[fk]:
                metafile = files[fk]['metadata_json']['path']
            if 'analytic' in files[fk]:
                image_file = files[fk]['analytic']['path']
            elif 'analytic_ntf' in files[fk]:
                image_file = files[fk]['analytic_ntf']['path']
            elif 'pansharpened' in files[fk]:
                image_file = files[fk]['pansharpened']['path']
            elif 'sr' in files[fk]:
                image_file = files[fk]['sr']['path']
            elif 'composite' in files[fk]:
                image_file = files[fk]['composite']['path']
                metafile = files[flist[-1]]['metadata']['path']
            meta = ac.planet.metadata_parse(metafile)
            if ('platform' in meta) & (os.path.exists(image_file)):
                input_type = 'Planet'
                break  ## exit loop
        except:
            pass ## continue to next sensor
        ## end Planet
        ################

        ################
        ## GF
        try:
            tiles, metafile = ac.gf.bundle_test(bundle)
            if metafile is not None:
                meta = ac.gf.metadata(metafile)
                if meta['SatelliteID'] in ['GF1', 'GF1D', 'GF6']:
                    input_type = 'GF'
                    break ## exit loop
        except:
            pass ## continue to next sensor
        ## end GF
        ################

        ################
        ## AMAZONIA
        try:
            files_xml, files_tiff = ac.amazonia.bundle_test(bundle)
            meta = ac.amazonia.metadata(files_xml[0])
            if meta['sensor'] in ['AMAZONIA1_WFI', 'CBERS4A_WFI']:
                input_type = 'AMAZONIA'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end AMAZONIA
        ################

        ################
        ## FORMOSAT
        try:
            tiles, metafile = ac.formosat.bundle_test(bundle)
            if metafile is not None:
                meta = ac.formosat.metadata(metafile)
                if meta['sensor'] in ['FORMOSAT5_RSI']:
                    input_type = 'FORMOSAT'
                    break ## exit loop
        except:
            pass ## continue to next sensor
        ## end FORMOSAT
        ################


        ################
        ## SDGSAT-1 KX10 MII
        try:
            metafiles, calfiles, imgfiles = ac.sdgsat.bundle_test(bundle)
            for mi, mf in enumerate(metafiles):
                meta = ac.sdgsat.metadata(mf)
                cal = ac.sdgsat.calibration(calfiles[mi])
                sensor = '{}_{}'.format(meta['SatelliteID'], meta['SensorID'])
                if sensor in ['KX10_MII']:
                    input_type = 'SDGSAT'
                    break ## exit loop
        except:
            pass ## continue to next sensor
        ## end SDGSAT-1 KX10 MII
        ################

        ################
        ## IKONOS2
        try:
            files_dict = ac.ikonos.bundle_test(bundle)
            metadata = ac.ikonos.metadata_parse(files_dict['metadata'])
            if 'Source Image Metadata' in metadata:
                if metadata['Source Image Metadata']['Sensor'][0] == 'IKONOS-2':
                    input_type = 'IKONOS'
                    break ## exit loop
            elif 'Product Order Metadata' in metadata:
                if metadata['Product Order Metadata']['Sensor Name'][0] == 'IKONOS-2':
                    input_type = 'IKONOS'
                    break ## exit loop
        except:
            pass ## continue to next sensor
        ## end IKONOS2
        ################


        ################
        ## DEIMOS2
        try:
            files_dict = ac.deimos.bundle_test(bundle)
            if 'MS' in files_dict:
                metadata = ac.deimos.metadata(files_dict['MS']['metadata'])
            elif 'PAN' in images:
                metadata = ac.deimos.metadata(files_dict['PAN']['metadata'])
            if metadata['MISSION'] == 'Deimos 2':
                    input_type = 'DEIMOS'
                    break ## exit loop
        except:
            pass ## continue to next sensor
        ## end DEIMOS2
        ################


        ################
        ## ECOSTRESS
        try:
            meta = ac.ecostress.attributes(bundle)
            if (meta['InstrumentShortName'] == 'ECOSTRESS') & (meta['ProcessingLevelID'] == 'L1B'):
                input_type = 'ECOSTRESS'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end ECOSTRESS
        ################

        ################
        ## ENMAP
        try:
            bundle_files = ac.enmap.bundle_test(bundle)
            metadata, band_data = ac.enmap.metadata_parse(bundle_files['METADATA'])
            if metadata['mission'].upper() + '_' + metadata['sensor'].upper() == 'ENMAP_HSI':
                input_type = 'ENMAP'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end ENMAP
        ################

        ################
        ## EMIT
        try:
            meta = ac.shared.nc_gatts(bundle)
            if (meta['instrument'] == 'EMIT') & ('L1B' in meta['title']):
                input_type = 'EMIT'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end EMIT
        ################

        ################
        ## SEVIRI
        try:
            ret = ac.seviri.bundle_test(bundle)
            meta = ac.seviri.metadata(ret)
            if (meta['AIID'] == 'SEVI'):
                input_type = 'SEVIRI'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end SEVIRI
        ################

        ################
        ## DIMAP
        try:
            dimfile, datfile = ac.dimap.bundle_test(bundle)
            meta = ac.dimap.metadata(dimfile)
            if (meta['met_format'] == 'DIMAP') & (meta['met_dataset'] == 'BEAM-PRODUCT'):
                input_type = 'DIMAP'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end DIMAP
        ################

        ################
        ## AVHRR
        try:
            image_file, meta_file = ac.avhrr.bundle_test(bundle)
            meta = ac.avhrr.metadata(meta_file)
            if ('AVHRR' in meta['sensor']):
                input_type = 'AVHRR'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end AVHRR
        ################

        ################
        ## HYPSO
        try:
            gatts = ac.shared.nc_gatts(bundle)
            datasets = ac.shared.nc_datasets(bundle)
            if gatts['instrument'].startswith('HYPSO-'):
                input_type = 'HYPSO'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end HYPSO
        ################

        ################
        ## Wyvern
        try:
            file, jf = ac.wyvern.bundle_test(bundle)
            meta, gatts = ac.wyvern.metadata_parse(jf)
            if gatts['sensor'].startswith('Wyvern'):
                input_type = 'Wyvern'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end Wyvern
        ################

        ################
        ## Haiyang
        try:
            ret = ac.haiyang.bundle_test(bundle)
            meta = ac.haiyang.metadata(ret['metadata'])
            if (meta['SatelliteID'] in ['HY-1C', 'HY-1D']):
                input_type = 'HAIYANG'
                break ## exit loop
        except:
            pass ## continue to next sensor
        ## end Haiyang
        ################

        ################
        break ## exit loop

    ## remove the extracted bundle if it could not be identified
    if (input_type is None) & (zipped) & (os.path.exists(bundle)) & (bundle != orig_bundle):
        shutil.rmtree(bundle)
        bundle = '{}'.format(orig_bundle)

    print('Identified {} as {} type'.format(orig_bundle, input_type))

    ## return input_type
    return(input_type, bundle, zipped, extracted_path)
