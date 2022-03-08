## def parse
## parses acolite settings files and merges with defaults
## written by Quinten Vanhellemont, RBINS
## 2021-03-09
## modifications:


def parse(sensor, settings=None, merge=True):
    import os, time
    import numpy as np
    import scipy.ndimage
    import acolite as ac

    ## read default settings for sensor
    if (sensor is not None) | (merge):
        setu = ac.acolite.settings.load(sensor)
    else:
        setu = {}

    ## add user settings
    if settings is not None:
        if type(settings) is str:
            sets = ac.acolite.settings.read(settings)
            if merge:
                for k in sets: setu[k] = sets[k]
            else:
                setu = {k:sets[k] for k in sets}
        elif type(settings) is dict:
            if merge:
                for k in settings: setu[k] = settings[k]
            else:
                setu = {k:settings[k] for k in settings}

    ## make sure luts setting is a list
    if 'luts' in setu:
        if type(setu['luts']) is not list: setu['luts'] = [setu['luts']]

    ## convert values from settings file into numbers
    int_list = ['s2_target_res', 'map_max_dim',
              'dsf_filter_box', 'dsf_tile_dimensions', 'dsf_intercept_pixels', 'dsf_smooth_box',
              'blackfill_wave', 'l2w_mask_wave', 'l2w_mask_cirrus_wave',
              'glint_force_band', 'glint_mask_rhos_wave', 'exp_wave1',  'exp_wave2',
              'l2w_mask_smooth_sigma',
              'flag_exponent_swir', 'flag_exponent_cirrus','flag_exponent_toa',
              'flag_exponent_negative', 'flag_exponent_outofscene',
              'rgb_red_wl','rgb_green_wl', 'rgb_blue_wl',
              'geometry_res', 'verbosity', 'map_dpi',
              'dsf_wave_range', 'l2w_mask_negative_wave_range', 'dsf_residual_glint_wave_range',
              'luts_pressures', 'nechad_range', 'dsf_minimum_segment_size',
              'netcdf_compression_level', 'netcdf_compression_least_significant_digit',
              'output_projection_xrange', 'output_projection_yrange']

    float_list = ['min_tgas_aot', 'min_tgas_rho',

                  'dsf_percentile', 'dsf_filter_percentile',
                  'dsf_min_tile_aot', 'dsf_min_tile_cover',

                  'exp_swir_threshold', 'exp_fixed_epsilon_percentile',
                  'exp_fixed_aerosol_reflectance_percentile',

                  'dem_pressure_percentile', 'elevation',
                  'uoz_default', 'uwv_default',
                  'wind_default', 'wind',

                  'blackfill_max', 'glint_mask_rhos_threshold',
                  'l2w_mask_threshold', 'l2w_mask_cirrus_threshold','l2w_mask_high_toa_threshold',
                  'map_auto_range_percentiles',
                  'rgb_min', 'rgb_max', 'rgb_autoscale_percentiles',
                  'gains_toa', 'offsets_toa',
                  'acstar3_max_wavelength', 'nechad_max_rhow_C_factor',
                  'output_projection_resolution', 'default_projection_resolution',
                  'worldview_reproject_resolution']

    ## convert values to numbers
    for k in setu:
        if k not in setu: continue
        if setu[k] is None: continue

        if type(setu[k]) is list:
            if k in int_list: setu[k] = [int(i) for i in setu[k]]
            if k in float_list: setu[k] = [float(i) for i in setu[k]]
        else:
            if k in int_list: setu[k] = int(setu[k])
            if k in float_list: setu[k] = float(setu[k])

    ## default pressure
    if 'pressure' in setu:
        setu['pressure'] = 1013.25 if setu['pressure'] is None else float(setu['pressure'])

    return(setu)
