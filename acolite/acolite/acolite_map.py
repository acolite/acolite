## def acolite_map
## function to create maps from NetCDF outputs
## written by Quinten Vanhellemont, RBINS
## 2021-03-11
## modifications: 2021-03-11 (QV) RGB outputs
##                2021-03-15 (QV) large update, including other parameters and mapping with pcolormesh
##                2021-04-01 (QV) changed plot_all option
##                2022-06-21 (QV) changed handling of l2_flags (if not int)
##                2022-11-14 (QV) changed rgb pcolormesh, added font options
##                2022-12-20 (QV) added different RGB stretches
##                2022-12-27 (QV) fix for RGB pcolormesh (data already stretched 0-1)
##                2022-12-29 (QV) fix for RGB raster (data scaled again to 255)
##                2022-12-30 (QV) moved rgb stretch to different function, added rgb gamma
##                                added scale bar option
##                2023-09-28 (QV) improved rgb mapping with mpl>3.7
##                2024-03-14 (QV) update settings handling
##                2024-04-16 (QV) use new gem NetCDF handling
##                2024-04-20 (QV) set out of range data to min/max for RGB scaling
##                2024-05-28 (QV) add gem open/close, generalised rgb outputs
##                2024-05-31 (QV) fix for generalised rgb outputs and the presence of other datasets containing rho
##                2024-11-20 (QV) changed RGB title labeling
##                2025-01-16 (QV) removed transformation to lower case for parameter names
##                                added user parameter scaling file
##                2025-02-04 (QV) improved settings handling, changed parameter label identification

def acolite_map(ncf, output = None,
                settings = None,
                plot_all = True,
                plot_datasets = [],
                plot_skip = ['lon', 'lat', 'l2_flags'],
                map_save = True,
                map_show = False ):

    import os, copy
    import numpy as np
    import acolite as ac

    import pyproj
    from osgeo import ogr, osr

    import matplotlib as mpl
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from PIL import Image

    ## get image x/y ratio
    #def get_imratio(im):
    #    imratio = im.shape[0]/im.shape[1]
    #    map_inches = (setu['map_max_dim']/setu['map_dpi'])
    #    if imratio > 1:
    #        figsize = map_inches, int(map_inches / imratio)
    #    else:
    #        figsize = int(map_inches / imratio), map_inches
    #    return(imratio, figsize)

    ## output map to file
    def output_map(im, par):
        rgb = len(im.shape) > 2
        norm, cmap = None, None

        ## find out parameter scaling to use
        if not rgb:
            ## get parameter config
            cparl = '{}'.format(par)
            sp = cparl.split('_')
            wave = None
            pard = None
            ## parameter is configured
            if cparl in pscale:
                pard = {k:pscale[cparl][k] for k in pscale[cparl]}
            else:
                ## is parameter configured with asterisk?
                for i in range(len(sp)):
                    park = '{}_*'.format('_'.join(sp[0:-i]))
                    if (park in pscale):
                        pard = {k:pscale[park][k] for k in pscale[park]}
                        wave = sp[-1]
                        break
            ## not configured
            if pard is None:
                pard = {'log':False, 'name':cpar, 'unit': ''}
                pard['cmap'] = setu['map_default_colormap']#'Planck_Parchment_RGB'
            ## do auto ranging
            if setu['map_auto_range'] | ('min' not in pard) | ('max' not in pard):
                drange = np.nanpercentile(im, setu['map_auto_range_percentiles'])
                pard['min'] = drange[0]
                pard['max'] = drange[1]

            ## parameter name and title
            part = '{}{} [{}]'.format(pard['name'], '' if wave is None else ' {} nm'.format(wave), pard['unit'])

            if pard['cmap'] == 'default': pard['cmap']=setu['map_default_colormap']
            ctfile = "{}/{}/{}.txt".format(ac.config['data_dir'], 'Shared/ColourTables', pard['cmap'])
            if os.path.exists(ctfile):
                pard['cmap'] = mpl.colors.ListedColormap(np.loadtxt(ctfile)/255.)
            else:
                try:
                    cmap = copy.copy(mpl.cm.get_cmap(pard['cmap']))
                except:
                    pard['cmap'] = setu['map_default_colormap']

            ## copy colour map to not set bad/under globally
            cmap = copy.copy(mpl.cm.get_cmap(pard['cmap']))
            cmap.set_bad(setu['map_fill_color'])
            if setu['map_fill_outrange']:
                cmap.set_under(setu['map_fill_color'])

            ## do log scaling
            if pard['log']:
                im = np.log10(im)
                pard['min'] = np.log10(pard['min'])
                pard['max'] = np.log10(pard['max'])
                part = 'log10 {}'.format(part)

            ## normalisation
            norm=mpl.colors.Normalize(vmin=pard['min'], vmax=pard['max'])#, clip=setu['map_fill_outrange'])
        else:
            part = r'$\rho_{}$ RGB'.format('{'+par.replace('rgb_rho', '')+'}')
            if setu['map_title_rgb_wavelengths']:
                part += ' ({})'.format(', '.join(['{:.0f} nm'.format(w) for w in rgb_used]))

        ## title and outputfile
        title = '{}\n{}'.format(title_base, part)
        ofile = '{}/{}_{}.{}'.format(odir, fn, par, setu['map_ext'])

        ## raster 1:1 pixel outputs
        if setu['map_raster']:
            if not rgb:
                ## scale to 255 int
                im = ac.shared.datascl(im, dmin=pard['min'], dmax=pard['max'])
                ## look up in colormap and convert to uint
                im = cmap.__call__(im)[:,:,0:3]*255
                im = im.astype(np.uint8)
            else:
                ## fix for new RGB scaling methods
                im *= 255
                im = im.astype(np.uint8)

            img = Image.fromarray(im)
            img.save(ofile)

        ## matplotlib outputs
        else:
            #imratio, figsize = get_imratio(im)
            figsize = None
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(1, 1, 1, projection=crs)
            plt.sca(ax)
            if crs is None:
                if setu['map_pcolormesh']:
                    if rgb: ## convert rgb values to color tuple before mapping
                        if mpl.__version__ > '3.7.0':
                            axim = plt.pcolormesh(lon, lat, im, shading='auto')
                        else:
                            mesh_rgb = im[:, :, :]
                            colorTuple = mesh_rgb.reshape((mesh_rgb.shape[0] * mesh_rgb.shape[1]), 3)
                            colorTuple = np.insert(colorTuple,3,1.0,axis=1)
                            axim = plt.pcolormesh(lon, lat, im[:,:,0], color=colorTuple, shading='auto')
                    else:
                        axim = plt.pcolormesh(lon, lat, im, norm=norm, cmap=cmap, shading='auto')
                        if scene_mask is not None:
                            plt.pcolormesh(lon, lat, scene_mask, cmap='gray', vmin=0, vmax=1)

                    plt.xlabel('Longitude (°E)', **font)
                    plt.ylabel('Latitude (°N)', **font)
                    if limit is not None:
                        plt.xlim(limit[1],limit[3])
                        plt.ylim(limit[0],limit[2])
                    plt.xticks(**font)
                    plt.yticks(**font)

                    ## rotate tick labels
                    plt.xticks(rotation = setu['map_xtick_rotation'])
                    plt.yticks(rotation = setu['map_ytick_rotation'])

                else:
                    axim = plt.imshow(im, norm=norm, cmap=cmap)
                    if scene_mask is not None:
                        plt.imshow(scene_mask,cmap='gray', vmin=0, vmax=1)

                    plt.axis('off')
            else: ## cartopy
                axim = ax.imshow(im, origin='upper', extent=img_extent, transform=image_crs, norm=norm, cmap=cmap)
                gl = ax.gridlines(draw_labels=True, linewidth=1, color=setu['map_gridline_color'], alpha=0.75, linestyle='--')
                gl.top_labels = False
                gl.left_labels = True
                gl.bottom_labels = True
                gl.right_labels = False
                ## set size and rotation of ticklabels
                gl.xlabel_style = {'size': setu['map_fontsize'], 'rotation': setu['map_xtick_rotation']}
                gl.ylabel_style = {'size': setu['map_fontsize'], 'rotation': setu['map_ytick_rotation']}

                ## format ticks
                from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER

            if setu['map_title']: plt.title(title, **font)

            ## add point markers
            if points is not None:
                for pname in points:
                    if 'px' not in points[pname]: continue
                    p = points[pname]
                    if 'facecolor' in p:
                        mfc = p['facecolor']
                        mew = 1.5
                    else:
                        mfc = None
                        mew = None
                    if 'edgecolor' in p:
                        mec = p['edgecolor']
                    else:
                        mec = None
                    if 'fontsize' in p:
                        fontsize = p['fontsize']
                    else:
                        fontsize = None
                    if 'markersize' in p:
                        markersize = p['markersize']
                    else:
                        markersize = None
                    ## plot marker
                    pplot = plt.plot(p['px'], p['py'], color=p['color'],
                                             marker=p['sym'], zorder=10,
                                             markersize = markersize, mec=mec, mfc=mfc, mew=mew)
                    ## plot marker label
                    if p['label_plot']:
                        text_color = 'white'
                        if mfc is not None: text_color = mfc
                        plt.text(p['pxl'], p['pyl'], p['label'], color=text_color, fontsize=fontsize, zorder=10,
                                 path_effects = [pe.withStroke(linewidth=2, foreground=p['color'])],
                                 verticalalignment=points[pname]['va'], horizontalalignment=points[pname]['ha'])
            ## end add point markers

            ## add the scalebar
            if setu['map_scalebar']:
                plt.plot(xsb, ysb, '-', linewidth = 2, color=setu['map_scalebar_color'], zorder=10)
                ## add the label
                plt.text(xsbl, ysbl, sclabel, color=setu['map_scalebar_color'], zorder=11,
                                 horizontalalignment='center', fontsize=setu['map_fontsize'])
            ## end scalebar

            ## color bars
            cbar = None
            if setu['map_colorbar']:
                if setu['map_colorbar_orientation'] == 'vertical':
                    if crs is None:
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes('right', size='5%', pad=0.05)
                    else:
                        cax = ax.inset_axes((1.02, 0, 0.02, 1)); #make a color bar axis
                    cbar = fig.colorbar(axim, cax=cax, orientation='vertical')
                    cbar.ax.set_ylabel(part)
                else:
                    if crs is None:
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes('bottom', size='5%', pad=0.05)
                    else:
                        cax = ax.inset_axes((0, -0.05, 1, 0.02)); #make a color bar axis
                    cbar = fig.colorbar(axim, cax=cax, orientation='horizontal')
                    cbar.ax.set_xlabel(part)

            #plt.tight_layout()

            ## we also make a colorbar for RGB, to keep map extents the same
            ## delete it here
            if rgb:
                if cbar is not None:
                    #cbar.solids.set_edgecolor("w")
                    #cbar.outline.set_visible(False)
                    #cbar.set_ticks([])
                    cbar.ax._visible = False

            if map_save:
                plt.savefig(ofile, dpi=setu['map_dpi'], bbox_inches='tight', facecolor='white')
                if setu['verbosity']>1: print('Wrote {}'.format(ofile))
            if map_show:
                plt.show()
            plt.close()

    ## open NetCDF
    opened = False
    if type(ncf) is str:
        gem = ac.gem.gem(ncf)
        opened = True
    else:
        gem = ncf

    ## get info from netcdf file
    datasets = [ds for ds in gem.datasets]
    imratio = None

    ## combine default and user defined settings
    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}
    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(gem.gatts['sensor'])
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults
    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    ## set font settings
    font = {'fontname':setu['map_fontname'], 'fontsize':setu['map_fontsize']}
    mpl.rc('text', usetex=setu['map_usetex'])

    scene_mask = None
    if ('l2_flags' in gem.datasets) & (setu['map_mask']):
        scene_mask = gem.data('l2_flags')
        ## convert scene mask to int if it is not (e.g. after reprojection)
        if scene_mask.dtype not in [np.int16, np.int32]: scene_mask = scene_mask.astype(int)
        scene_mask = (scene_mask & (2**setu['flag_exponent_outofscene']))
        scene_mask[scene_mask != (2**setu['flag_exponent_outofscene'])] = 0 #np.nan
        scene_mask[scene_mask == (2**setu['flag_exponent_outofscene'])] = 1
        scene_mask = scene_mask.astype(np.float32)
        scene_mask[scene_mask==0] = np.nan

    ## get output limit
    limit = None
    if setu['map_limit'] is not None:
        if type(setu['map_limit']) is list:
            limit = [float(l) for l in setu['map_limit']]
        else:
            limit = setu['limit']

    ## load parameter configuration
    pscale = ac.acolite.parameter_scaling(file = ac.config['parameter_labels'])
    pscale_u = ac.acolite.parameter_scaling(file = ac.config['parameter_labels_user'])
    for k in pscale_u: pscale[k] = pscale_u[k]

    crs = None
    if setu['map_projected']:
        if setu['map_cartopy']:
            try:
                import cartopy.crs as ccrs
                d, att = gem.data(gem.gatts['projection_key'], attributes = True)
                crs_proj = pyproj.Proj(att['crs_wkt'])
                img_extent = gem.gatts['xrange'][0],gem.gatts['xrange'][1],gem.gatts['yrange'][1],gem.gatts['yrange'][0]
                pcrs = pyproj.CRS(att['crs_wkt'])
                uzone = pcrs.coordinate_operation.name.split()[-1]
                zone = int(uzone[0:2])
                south = uzone[-1].upper() == 'S'
                crs = ccrs.UTM(zone, southern_hemisphere=south)
                image_crs = ccrs.UTM(zone, southern_hemisphere=south)
            except:
                print('Could not determine projection for cartopy')
                crs = None
                if ('lon' in gem.datasets) & ('lat' in gem.datasets): setu['map_pcolormesh'] = True
        else:
            crs = None
            if ('lon' in gem.datasets) & ('lat' in gem.datasets): setu['map_pcolormesh'] = True

    rhos_ds = [ds for ds in gem.datasets if 'rhos_' in ds]
    rhos_wv = [int(ds.split('_')[-1]) for ds in rhos_ds]

    bn = os.path.basename(ncf)
    fn = bn.replace('.nc', '')
    if 'output' in setu: output = setu['output']
    odir = os.path.dirname(ncf) if output is None else output
    if not os.path.exists(odir):
        os.makedirs(odir)

    isodate = ''
    if 'isodate' in gem.gatts: isodate = gem.gatts['isodate']
    elif 'time_coverage_end' in gem.gatts: isodate = gem.gatts['time_coverage_end']
    elif 'time_coverage_start' in gem.gatts: isodate = gem.gatts['time_coverage_start']

    if 'satellite_sensor' in gem.gatts:
        title_base = '{} {}'.format(gem.gatts['satellite_sensor'].replace('_', '/'), isodate.replace('T', ' ')[0:19])
    else:
        title_base = '{} {}'.format(gem.gatts['sensor'].replace('_', '/'), isodate.replace('T', ' ')[0:19])

    ## parameters to plot
    plot_parameters = []
    if plot_all:
        plot_parameters = [k for k in gem.datasets if k not in plot_skip]
    else:
        plot_parameters = [k for k in gem.datasets if k in plot_datasets]
    for k in setu:
        if k[0:7] != 'rgb_rho': continue
        if setu[k]: plot_parameters+=[k]

    ## handle wildcards
    for par in plot_parameters:
        if '*' in par:
            plot_parameters.remove(par)
            plot_parameters += [ds for ds in gem.datasets if ds[0:par.find('*')] == par[0:par.find('*')]]
    if len(plot_parameters) == 0: return

    ## load lat and lon
    if setu['map_pcolormesh'] | (setu['map_points'] is not None) | (setu['map_scalebar']):
        lon = gem.data('lon')
        lat = gem.data('lat')
        ## find region mid point
        mid = int(lat.shape[0]/2), int(lat.shape[1]/2)
        ## find lon and lat ranges (at mid point)
        lonw, lone = lon[mid[0], 0], lon[mid[0], -1]
        latn, lats = lat[0, mid[1]], lat[-1, mid[1]]
        lonr, latr = lone - lonw, latn - lats
        ## compute distance in one degree at mid point latitude
        lond, latd = ac.shared.distance_in_ll(lat=lat[mid[0], mid[1]])
        dd = lond*abs(lone-lonw)
        ## ratio of lon to lat distance
        lrat = lond/latd

    ## check if we need to add points to the map
    if 'map_points' in setu:
        if setu['map_points'] is None:
            points = None
        elif type(setu['map_points']) is dict:
            points = {p: setu['map_points'][p] for p in setu['map_points']}
        else:
            points = ac.shared.read_points(setu['map_points'])
        if type(points) is dict:
            ## find point and label positions
            for pname in points:
                p = points[pname]
                plon = float(p['lon'])
                plat = float(p['lat'])
                ## remove labeling tags if already present
                for k in ['px', 'py', 'ha', 'va', 'pxl', 'pyl']:
                    if k in points[pname]: del points[pname][k]
                if plon > np.nanmax(lon): continue
                if plon < np.nanmin(lon): continue
                if plat > np.nanmax(lat): continue
                if plat < np.nanmin(lat): continue
                if setu['map_pcolormesh'] & setu['map_projected']:
                    px, py = plon, plat
                elif setu['map_cartopy'] & setu['map_projected']:
                    px, py = crs_proj(plon, plat)
                else:
                    tmp = ((lon - plon)**2 + (lat - plat)**2)**0.5
                    i, j = np.where(tmp == np.nanmin(tmp))
                    py, px = i[0], j[0]
                ## track x and y
                points[pname]['px'] = px
                points[pname]['py'] = py
                ## get label position
                if p['label_side']=='right':
                    points[pname]['va']='center'
                    points[pname]['ha']='left'
                    xo, yo =  lonr * 0.02, latr * 0.0 #/ lrat

                if p['label_side']=='left':
                    points[pname]['va']='center'
                    points[pname]['ha']='right'
                    xo, yo =  lonr * -0.02, latr * 0.0 #/ lrat

                if p['label_side']=='top':
                    points[pname]['va']='bottom'
                    points[pname]['ha']='center'
                    xo, yo =  lonr * 0.0, latr * 0.02 #/ lrat

                if p['label_side']=='bottom':
                    points[pname]['va']='top'
                    points[pname]['ha']='center'
                    xo, yo =  lonr * 0.0, latr * -0.02 #/ lrat

                ## store label points
                if not setu['map_pcolormesh']:
                    tmp = ((lon - (plon+xo))**2 + (lat - (plat+yo))**2)**0.5
                    i, j = np.where(tmp == np.nanmin(tmp))
                    points[pname]['pxl'] = j[0]
                    points[pname]['pyl'] = i[0]
                else:
                    points[pname]['pxl'] = px + xo
                    points[pname]['pyl'] = py + yo
        else:
            points = None
    ## end points

    ## prepare scale bar
    ## approximate distance in one degree of longitude
    if setu['map_scalebar']:
        if setu['map_scalebar_position'] not in ['UR','UL','LL','LR']:
            print('Map scalebar position {} not recognised.')
            print('Using default map_scalebar_position=UL.'.format(scalepos))
            setu['map_scalebar_position'] = 'UL'
        posv = {'U': 0.85, 'L': 0.10}
        posh = {'R': 0.95, 'L': 0.05}
        if 'map_scalebar_position_vu' in setu:
            posv['U'] = float(setu['map_scalebar_position_vu'])
        if 'map_scalebar_position_vl' in setu:
            posv['L'] = float(setu['map_scalebar_position_vl'])
        if 'map_scalebar_position_hr' in setu:
            posh['R'] = float(setu['map_scalebar_position_hr'])
        if 'map_scalebar_position_hl' in setu:
            posh['L'] = float(setu['map_scalebar_position_hl'])

        if setu['map_scalebar_position'][0]=='U':
            latsc = lats+abs(latn-lats)*posv['U'] #0.87
        if setu['map_scalebar_position'][0]=='L':
            latsc = lats+abs(latn-lats)*posv['L'] #0.08
        if setu['map_scalebar_position'][1]=='R':
            lonsc = lonw+abs(lone-lonw)*posh['R'] #0.92
            scale_sign=-1.
        if setu['map_scalebar_position'][1]=='L':
            lonsc = lonw+abs(lone-lonw)*posh['L'] #0.08
            scale_sign=1.
        ## scale bar width
        if setu['map_scalebar_length'] is not None:
            scalelen = int(setu['map_scalebar_length'])
            unit = 'km'
        else:
            ## compute optimal scale length (as maximum fraction of image width)
            scalelen = dd * setu['map_scalebar_max_fraction']
            scalelen, unit = ac.shared.scale_dist(scalelen)
        ## compute scaleline
        sf = 1
        if unit == 'm':
            scaleline = (scalelen / 1000) / lond
        else:
            scaleline = scalelen / lond
        ## scalebar label
        sclabel = '{} {}'.format(scalelen*sf, unit)
        ## compute scalebar position
        xsb = (lonsc, lonsc+scale_sign*scaleline)
        ysb = (latsc,latsc)
        ## compute scalebar label position
        xsbl = lonsc + (scale_sign*scaleline)/2
        ysbl = latsc + (latr * 0.03)
        ## compute positions in pixels
        if not (setu['map_pcolormesh'] | setu['map_cartopy']):
            tmp = ((lon - xsb[0])**2 + (lat - ysb[0])**2)**0.5
            il, jl = np.where(tmp == np.nanmin(tmp))
            tmp = ((lon - xsb[1])**2 + (lat - ysb[1])**2)**0.5
            ir, jr = np.where(tmp == np.nanmin(tmp))
            ysb = (il[0], ir[0])
            xsb = (jl[0], jr[0])
            tmp = ((lon - xsbl)**2 + (lat - ysbl)**2)**0.5
            ip, jp = np.where(tmp == np.nanmin(tmp))
            xsbl, ysbl = jp, ip
        if setu['map_cartopy'] & setu['map_projected']:
            xsb, ysb = crs_proj(xsb, ysb)
            xsbl, ysbl = crs_proj(xsbl, ysbl)
    ## end prepare scale bar

    ## make plots
    for cpar in plot_parameters:
        if 'projection_key' in gem.gatts:
            if cpar in ['x', 'y', gem.gatts['projection_key']]: continue

        cparl = '{}'.format(cpar)
        ## RGB
        if (cpar[0:7] == 'rgb_rho'):
            ## find datasets for RGB compositing
            rgb_wave = [setu['rgb_red_wl'],setu['rgb_green_wl'],setu['rgb_blue_wl']]
            ds_base = [ds.split('_')[0:-1] for ds in gem.datasets if cpar.split('_')[1] in ds[0:len(cpar.split('_')[1])]]
            if len(ds_base) == 0:
                ds_base = '{}_'.format(cpar.split('_')[1])
            else:
                ds_base = '_'.join(ds_base[0]) + '_'
                if setu['add_band_name']: ds_base = '{}_'.format(cpar.split('_')[1])
            rho_ds = [ds for ds in gem.datasets if ds_base in ds[0:len(ds_base)]]
            rho_wv = [int(ds.split('_')[-1]) for ds in rho_ds]
            if len(rho_wv) < 2: continue

            ## read and stack rgb
            rgb_used = []
            for iw, w in enumerate(rgb_wave):
                wi, ww = ac.shared.closest_idx(rho_wv, w)
                rgb_used.append(ww)
                #ds_name = '{}{}'.format(ds_base,ww)
                ds_name = [ds for ds in gem.datasets if (ds_base in ds) & ('{:.0f}'.format(ww) in ds)][0]
                data = gem.data(ds_name, mask = False)

                ## autoscale rgb to percentiles
                if setu['rgb_autoscale']:
                    bsc = np.nanpercentile(data.data, setu['rgb_autoscale_percentiles'])
                else:
                    bsc = np.asarray((setu['rgb_min'][iw], setu['rgb_max'][iw]))
                gamma = setu['rgb_gamma'][iw]

                ## do RGB stretch
                tmp = ac.shared.rgb_stretch(data, gamma = gamma, bsc = bsc, stretch=setu['rgb_stretch'])

                ## mask
                #tmp[np.isnan(data)] = 1

                ## set min/max to rgb stretch
                tmp[data<=bsc[0]] = 0
                tmp[data>=bsc[1]] = 1

                ## stack RGB
                if iw == 0:
                    im = tmp
                else:
                    im = np.dstack((im, tmp))
        ## other parameters
        else:
            if cparl not in datasets:
                print('{} not in {}'.format(cpar, ncf))
                continue
            ds = [ds for di, ds in enumerate(gem.datasets) if cparl==datasets[di]][0]
            ## read data
            im = gem.data(ds)

        ## plot figure
        output_map(im, cpar)

    if opened:
        gem.close()
        del gem
