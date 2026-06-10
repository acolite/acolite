## def rgb_stack
## function to create RGB visualisation based on input array
## with fixed and autoscale options
##
## written by Quinten Vanhellemont, RBINS
## 2026-06-10
## modifications:

def rgb_stack(wavelength, data, scale = True,
              auto_range = False, auto_range_min = [5]*3, auto_range_max = [95]*3, auto_range_sigma = 0,
              rgb_wave = [650, 560, 480], rgb_min = [0] * 3, rgb_max = [0.15] * 3, ):
    import acolite as ac
    import numpy as np

    ## run through wavelengths
    for iw, w in enumerate(rgb_wave):
        wi, ws = ac.shared.closest_idx(wavelength, w)

        if scale:
            ## output range for plotting
            orange = (0,1)

            ## determine input range
            if not auto_range:
                irange = (rgb_min[iw], rgb_max[iw])
            else:
                if auto_range_sigma != 0:
                    mean = np.nanmean(data[wi, :,:])
                    std = np.nanstd(data[wi, :,:])
                    irange = (mean - std * auto_range_sigma, mean + std * auto_range_sigma)
                else:
                    irange = (np.nanpercentile(data[wi, :,:], auto_range_min[iw]),
                              np.nanpercentile(data[wi, :,:], auto_range_max[iw]))

            ## interpolate from irange to orange (0,1)
            d = np.interp(data[wi, :, :], irange, orange)
        else:
            ## no scaling
            d = data[wi, :, :] * 1.0

        ## stack data
        if iw == 0:
            im = d
        else:
            im = np.dstack((im, d))
    return(im)
