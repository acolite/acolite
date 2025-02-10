## def l1_convert
## converts BEAM/DIMAP file to l1r NetCDF for acolite
## currently only tested for S3 OLCI data
## written by Quinten Vanhellemont, RBINS
## 2023-02-14
## modifications:
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-16 (QV) use new gem NetCDF handling
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output=None, settings=None):
    import os, glob
    import numpy as np
    import scipy.interpolate
    import xml.etree.ElementTree as ET
    import time, datetime, dateutil.parser
    import acolite as ac

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
    if verbosity > 1: print('Starting conversion of {} scene{}'.format(nscenes, 's' if nscenes==1 else ''))

    ## start conversion
    ofile = None
    ofiles = []
    for bundle in inputfile:
        dimfile, datfile = ac.dimap.bundle_test(bundle)
        meta = ac.dimap.metadata(dimfile)

        t0 = time.time()
        limit = None
        sub = None
        vname = None

        dn = os.path.dirname(dimfile)
        bn = os.path.basename(dimfile)
        bn, ex = os.path.splitext(bn)

        sensor = None
        if bn[0:7] == 'S3A_OL_': sensor = 'S3A_OLCI'
        if bn[0:7] == 'S3B_OL_': sensor = 'S3B_OLCI'

        if sensor is None:
            print('Sensor not identified from file {}'.format(dimfile))
            print('DIMAP processing not configured for {}'.format(dimfile))
            continue

        rsrd = ac.shared.rsr_dict(sensor)
        rsr_bands = rsrd[sensor]['rsr_bands']

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        verbosity = setu['verbosity']
        if output is None: output = setu['output']
        if output is None: output = os.path.dirname(bundle)

        if setu['smile_correction']:
            print('Smile correction not implemented in DIMAP processing.')

        #data_dir = '{}/{}.data/'.format(dn, bn)
        img_bands = glob.glob('{}/*.img'.format(datfile))
        img_bands.sort()
        tpg_bands = glob.glob('{}/tie_point_grids/*.img'.format(datfile))
        tpg_bands.sort()

        ## tpg attributes from metadata
        ## to get tpg spacing and offset
        tpg_atts = meta['tpg_atts']

        ## parse date time
        dtime = dateutil.parser.parse(meta['start_date'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## read data
        data = {}
        data_tpg = {}
        for img in img_bands+tpg_bands:
            ib = os.path.basename(img)
            ib, ie = os.path.splitext(ib)

            ## read data
            cdata = ac.shared.read_band(img)

            ## read header data
            bhdr = '{}/{}.hdr'.format(os.path.dirname(img), ib)
            bhead = ac.shared.hdr(bhdr)

            ## convert to radiance
            cdata = cdata.astype(np.float32) * bhead['data gain values']
            cdata += bhead['data offset values']

            if img in img_bands:
                data[ib] = cdata
            if img in tpg_bands:
                data_tpg[ib] = cdata

        ## data shapes
        data_shape = data[list(data.keys())[0]].shape
        tpg_shape = data_tpg[list(data_tpg.keys())[0]].shape

        ## interpolate TPG
        ## 1d for interp2d
        subx = np.arange(0,data_shape[1])+0.5
        suby = np.arange(0,data_shape[0])+0.5
        ## 2d for RGI
        subx_ = np.tile(subx, (len(suby),1))
        suby_ = (np.tile(suby, len(subx)).reshape(subx_.shape[1], subx_.shape[0])).T
        tpg_int = {}
        tpg_nc = {}
        for tp in ['TP_latitude', 'TP_longitude', 'OAA', 'OZA', 'SAA', 'SZA',
                   'sea_level_pressure', 'total_columnar_water_vapour', 'total_ozone']:
            dtype_in = data_tpg[tp].dtype
            tpx = (np.arange(data_tpg[tp].shape[1])*float(tpg_atts[tp]['STEP_X'])) + float(tpg_atts[tp]['OFFSET_X'])
            tpy = (np.arange(data_tpg[tp].shape[0])*float(tpg_atts[tp]['STEP_Y'])) + float(tpg_atts[tp]['OFFSET_Y'])
            rgi = scipy.interpolate.RegularGridInterpolator((tpy, tpx), data_tpg[tp].astype(np.float64), bounds_error=False, fill_value=None)
            tpg_int[tp] = rgi((suby_,subx_)).astype(dtype_in)

        ## RAA
        tpg_int['RAA'] = np.abs(tpg_int['SAA']-tpg_int['OAA'])
        tpg_int['RAA'][tpg_int['RAA']>180]=np.abs(360-tpg_int['RAA'][tpg_int['RAA']>180])

        ## average geometry
        sza = np.nanmean(tpg_int['SZA'])
        vza = np.nanmean(tpg_int['OZA'])
        saa = np.nanmean(tpg_int['SAA'])
        vaa = np.nanmean(tpg_int['OAA'])
        raa = np.nanmean(tpg_int['RAA'])

        ## center lat and lon
        #clat = np.nanmean(tpg_int['TP_latitude'])
        #clon = np.nanmean(tpg_int['TP_longitude'])

        ## compute gas transmittance
        if setu['use_supplied_ancillary']:
            ## convert ozone from kg.m-2 to cm.atm
            setu['uoz_default'] = np.nanmean(tpg_int['total_ozone'])/0.02141419
            ## convert water from kg.m-2 to g.cm-2
            setu['uwv_default'] = np.nanmean(tpg_int['total_columnar_water_vapour'])/10
            setu['pressure'] = np.nanmean(tpg_int['sea_level_pressure'])
            print(setu['uoz_default'], setu['uwv_default'], setu['pressure'])

        gatts = {'sensor':sensor, 'sza':sza, 'vza':vza, 'saa':saa, 'vaa':vaa, 'raa': raa,
                     'isodate':isodate, 'global_dims':data_shape,
                     'se_distance': se_distance, 'acolite_file_type': 'L1R'}

        if limit is not None: gatts['limit'] = limit
        if sub is not None: gatts['sub'] = sub

        stime = dateutil.parser.parse(gatts['isodate'])
        oname = '{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## read rsr
        waves_mu = rsrd[sensor]['wave_mu']
        waves_names = rsrd[sensor]['wave_name']
        dnames = ['{}_radiance'.format(b) for b in rsr_bands]
        bnames = [b for b in rsr_bands]

        ## get F0 - not stricty necessary if using USGS reflectance
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsrd[sensor]['rsr'])

        ## add band info to gatts
        for bi, b in enumerate(rsr_bands):
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]
            gatts['{}_f0'.format(b)] = f0_b[b]

        ## write L1R netcdf
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = gatts

        ## output lat/lon
        if setu['output_geolocation']:
            if verbosity > 1: print('Writing geolocation')
            dsets = {'lon':'longitude', 'lat':'latitude'}
            for ds in dsets: gemo.write(ds, data[dsets[ds]])

        ## output geometry
        if setu['output_geometry']:
            if verbosity > 1: print('Writing geometry')
            dsets = {'sza':'SZA', 'vza':'OZA', 'raa': 'RAA',
                     'saa': 'SAA', 'vaa':'OAA', 'pressure': 'sea_level_pressure'}
            for ds in dsets: gemo.write(ds, tpg_int[dsets[ds]])

        ## is the F0 already corrected for Sun - Earth distance?
        ## see OLCI_L2_ATBD_Instrumental_Correction.pdf
        #se = ac.shared.distance_se(doy)
        se = 1.
        se2 = se**2

        ## cosine of solar zenith angle
        mus = np.cos(np.radians(tpg_int['SZA']))

        ## write rhot
        if verbosity > 1: print('Computing TOA reflectance')
        for iw, band in enumerate(rsr_bands):
            bidx = int(band[2:])
            print(band, bidx)
            wave = waves_names[band]
            ds = 'rhot_{}'.format(wave)
            if verbosity > 2: print('{} - Writing TOA data for {} nm'.format(datetime.datetime.now().isoformat()[0:19], wave), end='\n')

            ds_att  = {'wavelength':float(wave)}
            #for key in ttg: ds_att[key]=ttg[key][bnames[iw]]

            ## apply gains
            if setu['gains']:
                cg = 1.0
                if len(setu['gains_toa']) == len(rsr_bands):
                    cg = float(setu['gains_toa'][iw])
                if verbosity > 2: print('Applying gain {:.5f} for {}'.format(cg, ds))
                data['{}_radiance'.format(band)]*=cg

            ## write toa radiance
            if setu['output_lt']:
                gemo.write('Lt_{}'.format(wave), data['{}_radiance'.format(band)], ds_att = ds_att)
                if verbosity > 2: print('Converting bands: Wrote {} ({})'.format('Lt_{}'.format(wave), data['{}_radiance'.format(band)].shape))

            ## compute toa reflectance
            d = (np.pi * data['{}_radiance'.format(band)] * se2) / (data['solar_flux_band_{}'.format(bidx)] * mus)
            gemo.write(ds, d, ds_att = ds_att)
            if verbosity > 2: print('Converting bands: Wrote {} ({})'.format(ds, d.shape))
            d = None

        ## clear data
        data = None
        gemo.close()

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if limit is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)
    return(ofiles, setu)
