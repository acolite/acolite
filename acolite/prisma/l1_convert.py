## def l1_convert
## converts PRISMA HDF file to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-07-14

def l1_convert(inputfile, output=None, verbosity=0, vname = ''):
    import numpy as np
    import h5py, dateutil.parser, os
    import acolite as ac

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## get F0 for radiance -> reflectance computation
    f0 = ac.shared.f0_get()

    ofiles = []
    for file in inputfile:
        #f = h5py.File(file, mode='r')
        h5_gatts = ac.prisma.attributes(file)

        waves_vnir = h5_gatts['List_Cw_Vnir']
        bands_vnir = ['{:.0f}'.format(w) for w in waves_vnir]
        fwhm_vnir = h5_gatts['List_Fwhm_Vnir']
        n_vnir = len(waves_vnir)

        waves_swir = h5_gatts['List_Cw_Swir']
        bands_swir = ['{:.0f}'.format(w) for w in waves_swir]
        fwhm_swir = h5_gatts['List_Fwhm_Swir']
        n_swir = len(waves_swir)

        waves = [w for w in waves_vnir] + [w for w in waves_swir]
        fwhm = [f for f in fwhm_vnir] + [f for f in fwhm_swir]
        waves_names = ['{:.0f}'.format(w) for w in waves]
        instrument = ['vnir']*n_vnir + ['swir']*n_swir
        band_index = [i for i in range(n_vnir)] + [i for i in range(n_swir)]

        band_names_vnir = ['vnir_{}'.format(b) for b in range(0, n_vnir)]
        band_names_swir = ['swir_{}'.format(b) for b in range(0, n_swir)]

        rsr_vnir = {'vnir_{}'.format(b): ac.shared.gauss_response(waves_vnir[b], fwhm_vnir[b], step=0.1) for b in range(0, n_vnir)}
        rsr_swir = {'swir_{}'.format(b): ac.shared.gauss_response(waves_swir[b], fwhm_swir[b], step=0.1) for b in range(0, n_swir)}

        band_names = band_names_vnir + band_names_swir
        band_rsr = {}
        for b in rsr_vnir: band_rsr[b] = {'wave': rsr_vnir[b][0]/1000, 'response': rsr_vnir[b][1]}
        for b in rsr_swir: band_rsr[b] = {'wave': rsr_swir[b][0]/1000, 'response': rsr_swir[b][1]}

        idx = np.argsort(waves)
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], band_rsr)

        idx = np.argsort(waves)
        bands = {}
        for i in idx:
            cwave = waves[i]
            if cwave == 0: continue
            swave = '{:.0f}'.format(cwave)
            bands[swave]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                           'wave_name':waves_names[i],
                           'width': fwhm[i],
                           'i':i, 'index':band_index[i],
                           'rsr': band_rsr[band_names[i]],
                           'f0': f0d[band_names[i]],
                           'instrument':instrument[i],}

        gatts = {}

        isotime = h5_gatts['Product_StartTime']
        time = dateutil.parser.parse(isotime)

        doy = int(time.strftime('%j'))
        d = ac.shared.distance_se(doy)

        ## get geometry from l2 file if present
        l2file = file.replace('PRS_L1_STD_OFFL_', 'PRS_L2C_STD_')
        if os.path.exists(l2file):
            with h5py.File(l2file, mode='r') as f:
                vza = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Observing_Angle'][:]
                raa = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Rel_Azimuth_Angle'][:]
                sza = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Solar_Zenith_Angle'][:]
                gatts['vza'] = np.nanmean(np.abs(vza))
                gatts['raa'] = np.nanmean(np.abs(raa))
                gatts['sza'] = np.nanmean(np.abs(sza))
        else:
            print('PRISMA processing only supported when L2 geometry is present.')
            print('Please put {} in the same directory as {}'.format(os.path.basename(l2file), os.path.basename(file)))
            continue

        src = 'HCO'
        read_cube = True
        with h5py.File(file, mode='r') as f:
            lat = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Geolocation Fields']['Latitude_SWIR'][:]
            lon = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Geolocation Fields']['Longitude_SWIR'][:]
            ## read bands in spectral order
            if read_cube:
                vnir_data = h5_gatts['Offset_Vnir'] + \
                            f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['VNIR_Cube'][:]/h5_gatts['ScaleFactor_Vnir']
                swir_data = h5_gatts['Offset_Swir'] + \
                            f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['SWIR_Cube'][:]/h5_gatts['ScaleFactor_Swir']
            ## read LOS vectors
            x_ = f['KDP_AUX']['LOS_Vnir'][:,0]
            y_ = f['KDP_AUX']['LOS_Vnir'][:,1]
            z_ = f['KDP_AUX']['LOS_Vnir'][:,2]

        ## get vza/vaa
        #dtor = np.pi/180
        #vza = np.arctan2(y_,x_)/dtor
        #vaa = np.arctan2(z_,np.sqrt(x_**2+y_**2))/dtor
        #vza_ave = np.nanmean(np.abs(vza))
        #vaa_ave = np.nanmean(np.abs(vaa))
        sza_ave = h5_gatts['Sun_zenith_angle']
        saa_ave = h5_gatts['Sun_azimuth_angle']
        cossza = np.cos(sza_ave*(np.pi/180.))
        vza_ave = 0
        vaa_ave = 0

        if 'sza' not in gatts: gatts['sza'] = sza_ave
        if 'vza' not in gatts: gatts['vza'] = vza_ave
        if 'saa' not in gatts: gatts['saa'] = saa_ave
        if 'vaa' not in gatts: gatts['vaa'] = vaa_ave

        if 'raa' not in gatts:
            raa_ave = abs(gatts['saa'] - gatts['vaa'])
            while raa_ave >= 180: raa_ave = abs(raa_ave-360)
            gatts['raa'] = raa_ave

        mu0 = np.cos(gatts['sza']*(np.pi/180))
        muv = np.cos(gatts['vza']*(np.pi/180))

        if output is None:
            odir = os.path.dirname(file)
        else:
            odir = output

        obase  = '{}_{}_L1R'.format('PRISMA',  time.strftime('%Y_%m_%d_%H_%M_%S'))
        if not os.path.exists(odir): os.makedirs(odir)
        ofile = '{}/{}.nc'.format(odir, obase)

        gatts['obase'] = obase
        gatts['sensor'] = 'PRISMA'
        gatts['isodate'] = time.isoformat()

        gatts['band_waves'] = [bands[w]['wave'] for w in bands]
        gatts['band_widths'] = [bands[w]['width'] for w in bands]

        ac.output.nc_write(ofile, 'lat', np.flip(np.rot90(lat)), new=True, attributes=gatts)
        ac.output.nc_write(ofile, 'lon', np.flip(np.rot90(lon)))
        if os.path.exists(l2file):
            ac.output.nc_write(ofile, 'sza', np.flip(np.rot90(sza)))
            ac.output.nc_write(ofile, 'vza', np.flip(np.rot90(vza)))
            ac.output.nc_write(ofile, 'raa', np.flip(np.rot90(raa)))

        ## write TOA data
        for bi, b in enumerate(bands):
            wi = bands[b]['index']
            i = bands[b]['i']
            print('Reading rhot_{}'.format(bands[b]['wave_name']))

            if bands[b]['instrument'] == 'vnir':
                if read_cube:
                    cdata_radiance = vnir_data[:,wi,:]
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)
                else:
                    cdata_radiance = h5_gatts['Offset_Vnir'] + \
                            f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['VNIR_Cube'][:,i,:]/h5_gatts['ScaleFactor_Vnir']
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)

            if bands[b]['instrument'] == 'swir':
                if read_cube:
                    cdata_radiance = swir_data[:,wi,:]
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)
                else:
                    cdata_radiance = h5_gatts['Offset_Swir'] + \
                            f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['SWIR_Cube'][:,i,:]/h5_gatts['ScaleFactor_Swir']
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)

            ds_att = {k:bands[b][k] for k in bands[b] if k not in ['rsr']}

            if True:
                ## write toa radiance
                ac.output.nc_write(ofile, 'Lt_{}'.format(bands[b]['wave_name']), np.flip(np.rot90(cdata_radiance)),
                              dataset_attributes = ds_att)



            ## write toa reflectance
            ac.output.nc_write(ofile, 'rhot_{}'.format(bands[b]['wave_name']), np.flip(np.rot90(cdata)),
                              dataset_attributes = ds_att)
        ofiles.append(ofile)
    return(ofiles)
