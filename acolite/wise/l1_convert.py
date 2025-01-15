## def l1_convert
## converts WISE data to l1r NetCDF for acolite
## written by Raphael Mabit, ISMER-UQAR
## 2023-12-12
import math


def l1_convert(inputfile, output = None, settings = {}, verbosity = 0):
    import numpy as np
    import re
    import datetime, dateutil.parser, os, copy
    import acolite as ac

    ## parse sensor specific settings
    setu = ac.acolite.settings.parse('WISE_HSI', settings=settings)

    if output is None: output = setu['output']
    verbosity = setu['verbosity']
    poly = setu['polygon']
    limit = setu['limit']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## check if ROI polygon is given
    clip, clip_mask = False, None
    if poly is not None:
        if os.path.exists(poly):
            try:
                limit = ac.shared.polygon_limit(poly)
                print('Using limit from polygon envelope: {}'.format(limit))
                clip = True
            except:
                print('Failed to import polygon {}'.format(poly))

    ## get F0 for radiance -> reflectance computation
    # TODO check version of the coddington to use (2021 or 2023 ?)
    #  How the .txt file are produced ?  what is their unit ?
    #  The original files from https://lasp.colorado.edu/lisird/data/tsis1_hsrs_p1nm say W m-2 nm-1
    #  There is a factor of 1e3 between ACOLITE FO .txt and the original data, so the unit in ACOLITE are likely (W m-2 um-1)
    f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])

    ofiles = []

    for bundle in inputfile:
        print(bundle)

        imagename, pixhdrfile, pixfile, gluhdrfile, glufile, navlogfile = ac.wise.bundle_test(bundle)

        ## read header
        try:
            header = ac.shared.hdr(pixhdrfile)
        # When non utf-8 character in hdr file
        except UnicodeDecodeError as e:
            print(e)
            newlines = ac.wise.helper.rm_illegal_chr_envi_hdr(pixhdrfile)
            newheaderfile = pixhdrfile.replace('.hdr', '.copy.hdr')
            with open(newheaderfile, 'w') as f:
                f.writelines(newlines)
            header = ac.shared.hdr(newheaderfile)

        ## Parse the descrption as dictonary
        pairs = header['description']
        # remove elements from the list that contain dict split character ':'
        for elem in pairs[:]:
            if not any(char in elem for char in ':'):
                pairs.remove(elem)

        description = {key: value for key, value in (pair.strip(';').split(': ') for pair in pairs)}

        ## read Navcor_sum.log
        with open(navlogfile, 'r') as f:
            lines = f.readlines()

        for line in lines:
            line = line.strip()

            if line.startswith('Start'):
                line_strip = re.sub(r'(?<= )Line.*', '', line.replace('Start at:', '')).strip()
                navstime = line_strip

            if line.startswith('End'):
                line_strip = re.sub(r'(?<= )Line.*', '', line.replace('End   at:', '')).strip()
                navetime = line_strip

            if line.startswith('Average Value'):
                line_strip = line.replace('Average Value', '').replace('+ ground pixel size', '').replace('+',
                                                                                                          '').strip()
                # _logger.debug(line_strip.split(' '))
                values = [float(item) for item in line_strip.split(' ') if item.strip() != '']

        navlog = {
            'roll':values[0],
            'pitch':values[1],
            'heading':values[2],
            'distance':values[3],
            'height':values[4],
            'easting':values[5],
            'norhting':values[6]}

        ## set up projection
        warp_to, dct_prj, sub = None, None, None
        try:
            ## get projection from image
            dct = ac.wise.helper.parse_mapinfo(header)
        except:
            print('Could not determine image projection')
            dct = None

        ## find crop
        if (limit is not None) and (dct is not None):
            dct_sub = ac.shared.projection_sub(dct, limit)
            if dct_sub['out_lon']:
                if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                continue
            if dct_sub['out_lat']:
                if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                continue
            sub = dct_sub['sub']

        if dct is not None:
            if sub is None:
                dct_prj = {k:dct[k] for k in dct}
            else:
                dct_prj = {k:dct_sub[k] for k in dct_sub}

                ## updated 2022-03-28
                xyr = [min(dct_prj['xrange']),
                       min(dct_prj['yrange']),
                       max(dct_prj['xrange']),
                       max(dct_prj['yrange']),
                       dct_prj['proj4_string']]

                ## warp settings for read_band
                res_method = 'near'
                warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1], res_method)

        ## date and time
        stime = dateutil.parser.parse(navstime)
        etime = dateutil.parser.parse(str(stime.date())+'T'+navetime)
        otime = (etime-stime).seconds
        time = stime + datetime.timedelta(seconds=otime/2)
        doy = time.strftime('%j')
        se_distance = ac.shared.distance_se(doy)

        ## collect global attributes
        gatts = {}
        gatts['isodate'] = time.isoformat()
        gatts['sensor'] = description['instrument'].replace(" ", "_")
        #gatts['version'] = meta['version']
        gatts['doy'] = doy
        gatts['se_distance'] = se_distance
        # obase  = '{}_{}_L1R'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'))
        obase = '{}_{}_{}_L1R'.format(gatts['sensor'], imagename, time.strftime('%Y_%m_%d_%H_%M_%S'))
        gatts['obase'] = obase

        ## add projection info
        if dct_prj is not None:
            pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
            for k in pkeys:
                if k in dct_prj: gatts[k] = copy.copy(dct_prj[k])

            ## if we are clipping to a given polygon get the clip_mask here
            if clip:
                clip_mask = ac.shared.polygon_crop(dct_prj, poly, return_sub=False)
                clip_mask = clip_mask.astype(bool) == False

        ## add band info
        gatts['band_waves'] = header['wavelength']
        gatts['band_widths'] = header['fwhm']
        #gatts['band_names'] = header['band names']

        ## make rsr and bands dataset
        rsr = ac.shared.rsr_hyper(gatts['band_waves'], gatts['band_widths'], step=0.1)
        rsrd = ac.shared.rsr_dict(rsrd={gatts['sensor']:{'rsr':rsr}})
        band_rsr = rsrd[gatts['sensor']]['rsr']
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], band_rsr)

        ## make bands dataset
        bands = {}
        for bi, b in enumerate(band_rsr):
            cwave = rsrd[gatts['sensor']]['wave_nm'][b]
            swave = '{:.0f}'.format(cwave)
            bands[b]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                       'wave_name':'{:.0f}'.format(cwave),
                       'width': gatts['band_widths'][bi],
                       'rsr': band_rsr[b],'f0': f0d[b]}

        ## Sun geometry
        gatts['saa'] = float(header['sun azimuth'])
        gatts['sza'] = float(header['sun elevation'])

        ## Viewing geometry
        # TODO ACOLITE seem to use only the scene central viewing geometry. How much error is introduced by this assumption ?
        #  To compute the viewing geometry for every pixel in the scene, we could use the .glu and NAVCOR.LOG files.
        import transforms3d as t3d

        # test = t3d.taitbryan.euler2axangle((navlog['heading']+360)*(np.pi/180), navlog['roll']*(np.pi/180), navlog['pitch']*(np.pi/180))
        # vaa = math.atan(math.sqrt(pow(test[0][2], 2)+pow(test[0][1], 2))/test[0][0]) * 180/np.pi

        # For now, we take the scene central viewing geometry and assume a nadir viewing sensor
        gatts['vaa'] = 0
        gatts['vza'] = 0

        if 'raa' not in gatts:
            raa_ave = abs(gatts['saa'] - gatts['vaa'])
            while raa_ave >= 180: raa_ave = abs(raa_ave-360)
            gatts['raa'] = raa_ave

        mu0 = np.cos(gatts['sza']*(np.pi/180))
        muv = np.cos(gatts['vza']*(np.pi/180))

        imagefile = pixfile

        if output is None:
            odir = os.path.dirname(imagefile)
        else:
            odir = output
        if not os.path.exists(odir): os.makedirs(odir)
        ofile = '{}/{}.nc'.format(odir, obase)

        new = True
        if dct_prj is not None:
            print('Computing and writing lat/lon')
            ## offset half pixels to compute center pixel lat/lon
            # TODO check pixel reference (upper left or center ?) in dct_prj and helper.map_info.
            #   seems to have an inconsistency
            dct_prj['xrange'] = dct_prj['xrange'][0]+dct_prj['pixel_size'][0]/2, dct_prj['xrange'][1]-dct_prj['pixel_size'][0]/2
            dct_prj['yrange'] = dct_prj['yrange'][0]+dct_prj['pixel_size'][1]/2, dct_prj['yrange'][1]-dct_prj['pixel_size'][1]/2
            ## compute lat/lon
            lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel = False)
            print(lat.shape)
            ac.output.nc_write(ofile, 'lat', lat, new = new, attributes = gatts)
            lat = None
            ac.output.nc_write(ofile, 'lon', lon)
            lon = None
            new = False

        ## read data cube (faster)
        # WISE image are very large (17 ~ 60 G) cannot load them entirely in memory
        read_cube = False
        if read_cube:
            print('Reading WISE pix image')
            cube = ac.shared.read_band(imagefile, sub = sub, warp_to = warp_to).astype(np.float32)
            cube[cube == header['data ignore value']] = np.nan
            print(cube.shape)

        ## write TOA data
        for bi, b in enumerate(bands):
            print('Computing rhot_{} for {}'.format(bands[b]['wave_name'], gatts['obase']))
            ds_att = {k: bands[b][k] for k in bands[b] if k not in ['rsr']}

            ## read data
            if read_cube:
                cdata_radiance = 1.0 * cube[bi, :, :]
            else:
                cdata_radiance = ac.shared.read_band(imagefile, bi+1, sub=sub, warp_to = warp_to).astype(np.float32)

                # 'data ignore value' is empty in the header, looking in the data 0 is used as nodata
                #cdata_radiance[cdata_radiance == header['data ignore value']] = np.nan
                cdata_radiance[cdata_radiance == 0] = np.nan

            ## compute radiance
            # There is a 'radiance scale factor'
            # The unit says uW/cm2/nm/sr * 1000 (representing the packed data) so dividing by 1000 gives uW/cm2/nm/sr
            #print(description['units'])
            cdata_radiance = cdata_radiance.astype(np.float32) / float(header['radiance scale factor'])

            # There is no offset
            #cdata_radiance += header['data offset values'][bi]

            if (clip) & (clip_mask is not None): cdata_radiance[clip_mask] = np.nan

            output_lt = setu['output_lt']
            if output_lt:
                ## write toa radiance
                ac.output.nc_write(ofile, 'Lt_{}'.format(bands[b]['wave_name']), cdata_radiance,
                                            attributes = gatts, dataset_attributes = ds_att, new = new)
                new = False

            # compute reflectance
            # The WISE radiance data is in [uW cm-2 nm-1 sr-1]
            # Scaling ACOLITE FO file [W m2 um-1] to the original TSIS data [W m2 nm-1] is * 1e-3
            # Scaling [W m2 nm-1] to [uW cm-2 nm-1 sr-1] is 1e2
            cdata = cdata_radiance * (np.pi * gatts['se_distance'] * gatts['se_distance']) / ((bands[b]['f0']*1e-3*1e2) * mu0)
            cdata_radiance = None

            ac.output.nc_write(ofile, 'rhot_{}'.format(bands[b]['wave_name']), cdata,
                                            attributes = gatts, dataset_attributes = ds_att, new = new)
            cdata = None
            new = False
        cube = None

        ofiles.append(ofile)

    return(ofiles, setu)
