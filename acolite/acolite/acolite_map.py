## def acolite_map
## function to create maps from NetCDF outputs
## written by Quinten Vanhellemont, RBINS
## 2021-03-11
## modifications: 2021-03-11 (QV) RGB outputs
##

def acolite_map(ncf, output=None,
                dpi = 300, ext = 'png',
                map_save = True, map_show = False, map_title = True,
                map_projected = False, map_inches = 4,
                rgb_wave = [666, 560, 490],
                rgb_min = [0.0]*3, rgb_max = [0.15]*3,
                rgb_rhot=True, rgb_rhos=True):

    import os
    import numpy as np
    import acolite as ac
    import pyproj
    from osgeo import ogr, osr
    import cartopy.crs as ccrs

    def get_imratio(im):
        imratio = im.shape[0]/im.shape[1]
        if imratio > 1:
            figsize = map_inches, int(map_inches / imratio)
        else:
            figsize = int(map_inches / imratio), map_inches
        return(imratio, figsize)

    import matplotlib.pyplot as plt

    datasets = ac.shared.nc_datasets(ncf)
    gatts = ac.shared.nc_gatts(ncf)
    imratio = None

    crs = None
    if map_projected:
        try:
            proj4_string = gatts['proj4_string']
            p = pyproj.Proj(proj4_string)
            img_extent = gatts['xrange'][0],gatts['xrange'][1],gatts['yrange'][0],gatts['yrange'][1]
            pcrs = pyproj.CRS(gatts['proj4_string'])
            uzone = pcrs.coordinate_operation.name.split()[-1]
            zone = int(uzone[0:2])
            south = uzone[-1].upper() == 'S'
            crs = ccrs.UTM(zone, southern_hemisphere=south)
            image_crs = ccrs.UTM(zone, southern_hemisphere=south)
        except:
            print('Could not determine projection for cartopy')
            crs = None

    rhos_ds = [ds for ds in datasets if 'rhos_' in ds]
    rhos_wv = [int(ds.split('_')[1]) for ds in rhos_ds]

    bn = os.path.basename(ncf)
    odir = os.path.dirname(ncf) if output is None else output
    fn = bn.replace('.nc', '')

    title_base = '{} {}'.format(gatts['sensor'].replace('_', '/'), gatts['isodate'].replace('T', ' '))
    if (rgb_rhot) | (rgb_rhos):
        for ds_base in ['rhot_', 'rhos_']:
            if (ds_base == 'rhot_') & (not rgb_rhot): continue
            if (ds_base == 'rhos_') & (not rgb_rhos): continue

            rho_ds = [ds for ds in datasets if ds_base in ds]
            rho_wv = [int(ds.split('_')[1]) for ds in rho_ds]
            if len(rho_wv) < 3: continue

            ## read and stack rgb
            for iw, w in enumerate(rgb_wave):
                wi, ww = ac.shared.closest_idx(rho_wv, w)
                data = ac.shared.nc_data(ncf, '{}{}'.format(ds_base,ww))

                tmp = ac.shared.datascl(data.data, dmin=rgb_min[iw], dmax=rgb_max[iw])
                tmp[data.mask] = 255
                if iw == 0:
                    rgb = tmp
                else:
                    rgb = np.dstack((rgb, tmp))

            if imratio is None:  imratio, figsize = get_imratio(rgb)

            ## plot figure
            par = '{}rgb'.format(ds_base)
            ofile = '{}/{}_{}.{}'.format(odir, fn, par, ext)

            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(1, 1, 1, projection=crs)
            plt.sca(ax)
            if crs is None:
                plt.imshow(rgb)
                plt.axis('off')
            else:
                ax.imshow(rgb, origin='upper', extent=img_extent, transform=image_crs)
                gl = ax.gridlines(draw_labels=True)
                gl.xlabels_top = False
                gl.ylabels_left = True
                gl.xlabels_bottom = True
                gl.ylabels_right = False

            if map_title:
                plt.title('{} {}'.format(title_base, r'$\rho_{}$ RGB'.format(ds_base[-2])))

            plt.tight_layout()
            if map_save:
                plt.savefig(ofile, dpi=dpi)#, bbox_inches='tight')
            if map_show:
                plt.show()
            plt.close()
