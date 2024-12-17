## def acpar
## function to get 6SV like atmospheric correction parameters from libRadtran
##
## this can be achieved through a combination of simulations:
## 1) black surface at SUR and TOA to get the path radiance and downward transmittances (sun -> surface)
## 2) black surface at SUR with sun zenith angle = view zenith angle to get upward transmittance (surface -> sensor)
## 3) non-zero surface albedo to get spherical albedo of the atmosphere, by using the upward radiance
##
## written by Quinten Vanhellemont, RBINS
## 2024-10-30
## modifications: 2024-10-31 (QV) added kwargs, integrated in ACOLITE
##                2024-11-01 (QV) added wind keyword
##                2024-11-04 (QV) added quiet, run_delete and run_override, added phi offset

def acpar(sza, vza, raa, ## geometry - required
          wavelength_range = [300, 2500], ## wavelength range
          albedo = 0.2, ## surface albedo to get the spherical albedo of the atmosphere
          o3 = 0.3, h2o = 1.5, ## ozone and water vapour
          wind = None, pcl = 0,  sal = 34.3, ## settings for Cox & Munk, used when wind is not None
          parallel = True, ## use multiprocessing to run sims in parallel
          return_data = False, ## return simulation outputs
          quiet = True,
          run_base = None, ## user specified run name
          run_delete = True, ## delete run files
          run_override = True, ## override run files
          append_cfg =  [], ## list of lines to append at the end of the setting file
          **kwargs, ## keyword arguments to directly modify libRadtran settings
          ):
    import numpy as np
    import copy, os, time, datetime, string
    import acolite as ac

    ## cosine view and sun angles
    cos_vza = np.cos(np.radians(vza))
    cos_sza = np.cos(np.radians(sza))

    ## set up run name
    if run_base is None:
        run_base =  datetime.datetime.now().strftime('%Y%m%d_%H%M%S') + \
                    '_{}'.format(''.join(np.random.choice(list(string.ascii_uppercase), size=6, replace=True)))

    ## run 1 - basic setup
    lr1 = ac.rtm.libradtran.lira(run_name = run_base + '_run1')
    lr1.override = run_override
    lr1.delete = run_delete
    lr1.quiet = quiet
    lr1.cfg['wavelength'] = '{} {}'.format(wavelength_range[0],wavelength_range[1])
    lr1.cfg['sza'] = sza
    lr1.cfg['phi0'] = '0'
    lr1.cfg['umu'] = cos_vza
    lr1.cfg['phi'] = 180 - raa ## for libradtran: 0 is looking towards the sun, 180 is sun in the back of the sensor
    lr1.cfg["mol_modify O3"] = '{} DU'.format(o3*1000)
    lr1.cfg["mol_modify H2O"] = '{} MM'.format(h2o*10)
    #lr1.cfg["mol_abs_param reptran"] = 'coarse'

    ## direct modification of some settings through kwargs
    ## some libRadtran setting have spaces and are here simplified
    ## (more robust method may be needed if finer control is wanted)
    for k in kwargs:
        if k == 'O3':
            lr1.cfg["mol_modify O3"] = '{} DU'.format(kwargs[k])
        elif k == 'H2O':
            lr1.cfg["mol_modify H2O"] = '{} MM'.format(kwargs[k])
        elif k == 'reptran':
            lr1.cfg["mol_abs_param reptran"] = '{}'.format(kwargs[k])
        elif (k == 'wavelength') & (type(kwargs[k]) == list):
            lr1.cfg[k] = '{} {}'.format(kwargs[k][0], kwargs[k][-1])
        elif k == 'vza':
            lr1.cfg["umu"] = '{}'.format(np.cos(np.radians(kwargs[k])))
        elif k == 'raa':
            lr1.cfg["phi"] = '{}'.format(kwargs[k])
        else:
            lr1.cfg[k] = '{}'.format(kwargs[k])
    ## append settings lines directly
    for l in append_cfg: lr1.cfg[l] = ''

    if not quiet:
        print('libRadtran settings: ')
        for c in lr1.cfg:
            print('\t{} {}'.format(c, lr1.cfg[c]))

    ## copy setup for next runs
    lr2 = copy.deepcopy(lr1)
    lr3 = copy.deepcopy(lr1)

    ## run 2 - reciprocal sza = vza
    lr2.run_name = run_base + '_run2'
    lr2.set_file_paths()
    lr2.cfg['sza'] = vza

    ## run 3 - run with non-zero albedo
    lr3.run_name = run_base + '_run3'
    lr3.set_file_paths()
    lr3.cfg['albedo'] = albedo

    processes = 3
    ## run 4, only if wind is set
    if wind is not None:
        u10 = np.max((1, np.float32(wind))) ## set to at least 1 m/s
        lr4 = copy.deepcopy(lr1)
        lr4.run_name = run_base + '_run4'
        lr4.set_file_paths()
        lr4.cfg['albedo'] = albedo
        lr4.cfg["brdf_cam u10"] = "{}".format(u10)
        lr4.cfg["brdf_cam pcl"] = "{}".format(pcl)
        lr4.cfg["brdf_cam sal"] = "{}".format(sal)
        processes = 4

    if not quiet: print('Running simulations')
    t0 = time.time()
    if parallel: ## mp runs
        from multiprocessing import Pool
        with Pool(processes = processes) as pool:
            if processes == 3:
                [lr1, lr2, lr3] = pool.map(ac.rtm.libradtran.execute, [lr1, lr2, lr3])
            if processes == 4:
                [lr1, lr2, lr3, lr4] = pool.map(ac.rtm.libradtran.execute, [lr1, lr2, lr3, lr4])
    else: ## sequential runs
        lr1.run()
        lr2.run()
        lr3.run()
        if processes == 4: lr4.run()

    if not quiet: print('Total time {:.1f}s'.format(time.time()-t0))

    ## extract results and delete classes
    d1 = lr1.result
    del lr1
    d2 = lr2.result
    del lr2
    d3 = lr3.result
    del lr3
    if processes == 4:
        d4 = lr4.result
        del lr4

    ## output dict
    out = {}
    if return_data:
        out['d1'] = d1
        out['d2'] = d2
        out['d3'] = d3

    ## compute parameters
    ## wavelength
    out['wavelength'] = d1['TOA']['lambda'] * 1.0

    ## path reflectance
    if ('uu' in d1['TOA']) & ('edir' in d1['TOA']):
        #out['rho_path'] = (d1['TOA']['uu'] * np.pi)/(d1['TOA']['edir'] * cos_sza)
        out['rho_path'] = (d1['TOA']['uu'] * np.pi)/(d1['TOA']['edir'])

    ## sun-surface transmittance
    out['Td_tot'] = d1['SUR']['eglo']/(d1['TOA']['eglo'])
    out['Td_dir'] = d1['SUR']['edir']/(d1['TOA']['eglo'])
    out['Td_dif'] = d1['SUR']['edn']/(d1['TOA']['eglo'])

    ## surface-sensor transmittance
    out['Tu_tot'] = d2['SUR']['eglo']/(d2['TOA']['eglo'])
    out['Tu_dir'] = d2['SUR']['edir']/(d2['TOA']['eglo'])
    out['Tu_dif'] = d2['SUR']['edn']/(d2['TOA']['eglo'])

    ## spherical albedo - should be equivalent
    ## based on surface irradiance for two surface albedos, e.g. 0 for run 1
    #sa = (d1['SUR']['eglo'] - d3['SUR']['eglo']) / (0.0 * d1['SUR']['eglo'] - albedo * d3['SUR']['eglo'])
    #sa = (d1['SUR']['eglo'] - d3['SUR']['eglo']) / (- albedo * d3['SUR']['eglo'])
    ## based on upward radiance just above the surface for non-zero albedo
    if ('uu' in d3['SUR']) & ('eglo' in d1['SUR']):
        out['sa'] = (d3['SUR']['uu'] - d1['SUR']['eglo']*albedo/np.pi) / (d3['SUR']['uu']*albedo)

    del d1
    del d2
    del d3

    ## add Cox & Munk
    if processes == 4:
        if return_data: out['d4'] = d4
        if ('uu' in d4['TOA']) & ('edir' in d4['TOA']):
            #out['rho_path_ocean'] = (d4['TOA']['uu']* np.pi)/(d4['TOA']['edir'] * cos_sza)
            out['rho_path_ocean'] = (d4['TOA']['uu']* np.pi)/(d4['TOA']['edir'])

        del d4

    return(out)
