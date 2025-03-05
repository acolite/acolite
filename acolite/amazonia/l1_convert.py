## def l1_convert
## converts AMAZONIA bundle to l1r NetCDF for acolite-gen
## note: average geometry from both cameras is used
##       average RSR from both cameras is used (1-2% difference in Blue, less in other bands)
##       view azimuth is unknown, taken as 0 degrees (probably OK for low vza data)
## written by Quinten Vanhellemont, RBINS
## 2022-01-14
## modifications: 2022-02-06 (QV) added vza
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-17 (QV) use new gem NetCDF handling
##                2025-01-30 (QV) moved polygon limit and limit buffer extension
##                2025-02-02 (QV) removed percentiles
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output = None, settings = None):
    import os, zipfile, shutil
    import dateutil.parser, time
    import numpy as np
    import acolite as ac
    import re

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    verbosity = setu['verbosity']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## start with last file in time
    inputfile.sort()
    inputfile.reverse()

    new = True
    warp_to = None

    ofile = None
    ofiles = []

    for bundle in inputfile:
        t0 = time.time()
        files_xml, files_tiff = ac.amazonia.bundle_test(bundle)

        ## read meta data
        meta = {}
        bands = []
        for fm in files_xml:
            bbn = os.path.splitext(os.path.basename(fm))[0]
            band = bbn[bbn.find('BAND')+4:]
            bands.append(band)
            meta[band] = ac.amazonia.metadata(fm)

        meta['sensor'] = meta[bands[0]]['sensor']
        dtime = dateutil.parser.parse(meta[bands[0]]['isotime'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## merge sensor specific settings
        if new:
            ## get sensor specific defaults
            setd = ac.acolite.settings.parse(meta['sensor'])
            ## set sensor default if user has not specified the setting
            for k in setd:
                if k not in ac.settings['user']: setu[k] = setd[k]
            ## end set sensor specific defaults

            verbosity = setu['verbosity']
            if output is None: output = setu['output']
            if output is None: output = os.path.dirname(bundle)

            extend_region = setu['extend_region']
            ## check if merging settings make sense
            if (setu['limit'] is None) & (setu['merge_tiles']):
                if verbosity > 0: print("Merging tiles not supported without ROI limit")
                setu['merge_tiles'] = False
            if setu['merge_tiles']: extend_region = True

        sub = None

        ## read rsr
        rsrf = ac.path+'/data/RSR/{}.txt'.format(meta['sensor'])
        rsr, rsr_bands = ac.shared.rsr_read(rsrf)
        waves = np.arange(250, 2500)/1000
        waves_mu = ac.shared.rsr_convolute_dict(waves, waves, rsr)
        waves_names = {'{}'.format(b):'{:.0f}'.format(waves_mu[b]*1000) for b in waves_mu}

        ## get F0 - not stricty necessary if using USGS reflectance
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsr)

        ## gains
        gains_dict = None
        if (setu['gains']) & (setu['gains_toa'] is not None):
            if len(setu['gains_toa']) == len(rsr_bands):
                gains_dict = {b: float(setu['gains_toa'][ib]) for ib, b in enumerate(rsr_bands)}

        saa = []
        for band in bands:
            if 'main' in meta[band]:
                saa.append(meta[band]['main']['saa'])
            else:
                saa.append(meta[band]['leftCamera']['saa'])
                saa.append(meta[band]['rightCamera']['saa'])
        sza = []
        for band in bands:
            if 'main' in meta[band]:
                sza.append(meta[band]['main']['sza'])
            else:
                sza.append(meta[band]['leftCamera']['sza'])
                sza.append(meta[band]['rightCamera']['sza'])

        if 'main' in meta[band]:
            vza = meta[band]['main']['vza']
        else:
            vza = meta[band]['leftCamera']['vza']/2 + meta[band]['rightCamera']['vza']/2
        vaa = 0

        gatts = {'sensor':meta['sensor'],
                 'satellite_sensor':meta['sensor'],
                 'isodate':isodate,
                 'sza': np.nanmean(np.asarray(sza)),
                 'saa': np.nanmean(np.asarray(saa)),
                 'vza':vza, 'vaa': vaa,
                 'se_distance': se_distance,
                 'acolite_file_type': 'L1R'}

        gatts['mus'] = np.cos(gatts['sza']*(np.pi/180.))
        if 'raa' not in gatts:
            raa_ave = abs(gatts['saa'] - gatts['vaa'])
            while raa_ave >= 180: raa_ave = abs(raa_ave-360)

            #asub = np.where(raa_ave >= 180)
            #raa_ave[asub] = abs(raa_ave[asub]-360)
            gatts['raa'] = raa_ave

        stime = dateutil.parser.parse(gatts['isodate'])
        oname = '{}_{}{}'.format(gatts['satellite_sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'), '_merged' if setu['merge_tiles'] else '')
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])

        ## output file information
        if (setu['merge_tiles'] is False) | (ofile is None):
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile
        elif (setu['merge_tiles']) & (ofile is None):
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile

        ## add band info to gatts
        for b in rsr_bands:
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]
            gatts['{}_f0'.format(b)] = f0_b[b]

        dct = ac.shared.projection_read(files_tiff[0])
        gatts['scene_xrange'] = dct['xrange']
        gatts['scene_yrange'] = dct['yrange']
        gatts['scene_proj4_string'] = dct['proj4_string']
        gatts['scene_pixel_size'] = dct['pixel_size']
        gatts['scene_dims'] = dct['dimensions']
        if 'zone' in dct: gatts['scene_zone'] = dct['zone']

        ## check crop
        if (sub is None) & (setu['limit'] is not None):
            dct_sub = ac.shared.projection_sub(dct, setu['limit'], four_corners=True)
            if dct_sub['out_lon']:
                if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                continue
            if dct_sub['out_lat']:
                if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                continue
            sub = dct_sub['sub']
        else:
            if extend_region:
                print("Can't extend region if no ROI limits given")
                extend_region = False

        ##
        if sub is None:
            dct_prj = {k:dct[k] for k in dct}
        else:
            gatts['sub'] = sub
            gatts['limit'] = setu['limit']
            ## get the target NetCDF dimensions and dataset offset
            if (warp_to is None):
                if (extend_region): ## include part of the roi not covered by the scene
                    dct_prj = {k:dct_sub['region'][k] for k in dct_sub['region']}
                else: ## just include roi that is covered by the scene
                    dct_prj = {k:dct_sub[k] for k in dct_sub}
        ## end cropped

        ## get projection info for netcdf
        if setu['netcdf_projection']:
            nc_projection = ac.shared.projection_netcdf(dct_prj, add_half_pixel=True)
        else:
            nc_projection = None

        ## save projection keys in gatts
        pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
        for k in pkeys:
            if k in dct_prj: gatts[k] = dct_prj[k]

        ## warp settings for read_band
        ## updated 2021-10-28
        xyr = [min(dct_prj['xrange']),
               min(dct_prj['yrange']),
               max(dct_prj['xrange']),
               max(dct_prj['yrange']),
               dct_prj['proj4_string']]

        res_method = 'average'
        warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1], res_method)

        ## store scene and output dimensions
        gatts['scene_dims'] = dct['ydim'], dct['xdim']
        gatts['global_dims'] = dct_prj['dimensions']

        ## new file for every bundle if not merging
        new = True
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.nc_projection = nc_projection
        new = False

        ## if we are clipping to a given polygon get the clip_mask here
        if setu['polygon_clip']:
            clip_mask = ac.shared.polygon_crop(dct_prj, setu['polygon'], return_sub=False)
            clip_mask = clip_mask.astype(bool) == False

        if (os.path.exists(ofile) & (not new)):
            datasets = gemo.datasets
        else:
            datasets = []

        ## write lat/lon
        if (setu['output_geolocation']):
            if ('lat' not in datasets) or ('lon' not in datasets):
                if verbosity > 1: print('Writing geolocation lon/lat')
                lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel=True)
                gemo.write('lon', lon)
                lon = None
                if verbosity > 1: print('Wrote lon')
                gemo.write('lat', lat)
                lat = None
                if verbosity > 1: print('Wrote lat')

        ## write x/y
        if (setu['output_xy']):
            if ('x' not in datasets) or ('y' not in datasets):
                if verbosity > 1: print('Writing geolocation x/y')
                x, y = ac.shared.projection_geo(dct_prj, xy=True, add_half_pixel=True)
                gemo.write('xm', x)
                x = None
                if verbosity > 1: print('Wrote x')
                gemo.write('ym', y)
                y = None
                if verbosity > 1: print('Wrote y')

        ## convert bands
        b_vaa = []
        b_vza = []
        for bi, b in enumerate(rsr_bands):
            bfile = files_tiff[bi]
            if 'BAND{}'.format(b) not in bfile:
                print('Cannot fine file for BAND{}'.format(b))
                continue

            ## read data
            md, data = ac.shared.read_band(bfile, warp_to=warp_to, gdal_meta=True)
            nodata = data == np.uint16(0)

            ## convert to Lt
            if 'main' in meta[b]:
                #data = data.astype(np.float32) * meta[b]['main']['{}_absoluteCalibrationCoefficient'.format(b)]
                data = data.astype(np.float32) * 0.4

                if False:
                    ## table 8 from Pinto et al 2016 doi:10.3390/rs8050405
                    dn_gain = {'13':0.44, '14':0.47, '15':0.37, '16':0.34}
                    dn_offset = {'13':-19, '14':8, '15':-4, '16':3}
                    ## 0 offset
                    dn_gain = {'13':0.379, '14':0.498, '15':0.360, '16':0.351}
                    dn_offset = {'13':0, '14':0, '15':0, '16':0}
                    data = (data.astype(np.float32) * dn_gain[b]) + dn_offset[b]
            else:
                data = data.astype(np.float32) * meta[b]['leftCamera']['{}_absoluteCalibrationCoefficient'.format(b)]

            ## convert to rhot
            f0 = gatts['{}_f0'.format(b)]/10
            data *= (np.pi * gatts['se_distance']**2) / (f0 * np.nanmean(gatts['mus']))

            data[nodata] = np.nan
            if (setu['polygon_clip']): data[clip_mask] = np.nan

            ds = 'rhot_{}'.format(waves_names[b])
            ds_att = {'wavelength':waves_mu[b]*1000}

            ## get angles from tiff metadata
            for k in ['Elevation', 'Azimuth', 'Gain', 'Integration Time']:
                match = re.search(r'{} (\d+(?:\.\d+)?)'.format(k), md['TIFFTAG_IMAGEDESCRIPTION'])
                if match: ds_att[k] = float(match.group(1))
            ds_att['Zenith'] = 90 - ds_att['Elevation']
            b_vaa.append(ds_att['Azimuth'])
            b_vza.append(ds_att['Zenith'])

            if setu['gains'] & (gains_dict is not None):
                ds_att['toa_gain'] = gains_dict[b]
                data *= ds_att['toa_gain']
                if verbosity > 1: print('Converting bands: Applied TOA gain {} to {}'.format(ds_att['toa_gain'], ds))

            ## write to netcdf file
            gemo.write(ds, data, replace_nan = True, ds_att = ds_att)
            if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))

        ## update geometry from band tags
        gatts['vza'] = np.nanmean(b_vza)
        gatts['vaa'] = np.nanmean(b_vaa)
        raa_ave = abs(gatts['saa'] - gatts['vaa'])
        while raa_ave >= 180: raa_ave = abs(raa_ave-360)
        gatts['raa'] = raa_ave
        ## update nc gatts
        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.gatts_update()
        gemo.close()

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if setu['limit'] is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
