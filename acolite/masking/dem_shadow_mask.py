## def dem_shadow_mask
## makes shadow mask for given dem and sun azimuth and zenith angles, as well as dem scale
## based on the Richens 1997 shadow volume algorithm "Image processing for urban scale environmental modelling, Paul Richens, 1997b"
## written by Quinten Vanhellemont, RBINS
## 2024-01-25
## modifications:
##

def dem_shadow_mask(dem, saa, sza, dem_scale, loop = False):
    import acolite as ac
    import numpy as np
    import dateutil.parser

    ## set up sun angles
    ## note saa should be adjusted for grid convergence
    azi = np.radians(np.abs(270-saa)) ## rotate axis to have north up - numpy defaults to row major?
    ele = np.radians(90-sza)
    pixel_scale = 1/dem_scale
    sizey, sizex = dem.shape
    max_dem = np.nanmax(dem)
    sin_azi, cos_azi, tan_azi = np.sin(azi), np.cos(azi), np.tan(azi)
    sign_sin_azi = np.trunc(np.sign(sin_azi))
    sign_cos_azi = np.trunc(np.sign(cos_azi))
    dssin = np.abs(1/sin_azi)
    dscos = np.abs(1/cos_azi)
    tan_ele_scale = np.tan(ele) / pixel_scale

    ## initial position
    index = 1
    dx = 0
    dy = 0
    dz = 0.0
    dz = 0

    ## array copy
    shdem = 1.0 * dem

    ## run through shadow volume algorithm
    while (dz < max_dem) & (np.abs(dx) < sizex) & (np.abs(dy) < sizey):
        if ((np.pi / 4 <= azi) & (azi < 3 * np.pi / 4)) | ((5 * np.pi / 4 <= azi) & (azi < 7 * np.pi / 4)):
            dy = sign_sin_azi * index
            dx = -1 * sign_cos_azi * np.abs(np.round(index / tan_azi))
            ds = 1 * dssin
        else:
            dy = sign_sin_azi * np.abs(np.round(index * tan_azi))
            dx = -1 * sign_cos_azi * index
            ds = 1 * dscos

        dz = ds * index * tan_ele_scale

        absdx = np.abs(dx)
        absdy = np.abs(dy)

        xc1 = int(((dx + absdx) / 2) + 1)
        yc1 = int(((dy + absdy) / 2) + 1)

        xp1 = int(-((dx - absdx) / 2) + 1)
        xp2 = int((sizex - (dx + absdx) / 2))

        yp1 = int(-((dy - absdy) / 2) + 1)
        yp2 = int((sizey - (dy + absdy) / 2))

        xc_offset = xc1 - xp1
        yc_offset = yc1 - yp1

        if loop: ## very slow
            for i in range(yp1, yp2):
                for j in range(xp1, xp2):
                    shdem[i, j] = np.max((shdem[i,j], dem[i+yc_offset, j+xc_offset]-dz))
        else: ## fast
            shdem[yp1:yp2, xp1:xp2] = np.nanmax((shdem[yp1:yp2, xp1:xp2],
                                                 shdem[yp1+yc_offset:yp2+yc_offset, xp1+xc_offset:xp2+xc_offset]-dz), axis=0)

        #if verbosity > 2: print(index, dx, dy, dz)
        index+=1

    ## compute shadow mask
    mask_shadow = 1.0*(shdem > dem)
    return(mask_shadow)
