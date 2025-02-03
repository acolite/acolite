## class ilut
## class to import ACOLITE luts and perform interpolation
## written by Quinten Vanhellemont, RBINS
## 2025-02-03
## modifications:

import numpy as np
import scipy.interpolate
import acolite as ac

class ilut(object):
    def __init__(self, lut_base = None, lut_models = None, lut_pressures = None, lut_base_interface = None,
                  lut_dir = None, pressure_range = None, get_remote = True, keep_lut = False,
                  lut_par_defaults = ['wl',  'utotr', 'dtotr', 'phar', 'asray', 'tray', 'rorayl',
                                    'utota', 'dtota', 'phaa', 'asaer', 'taer', 'roaero',
                                    'utott', 'dtott', 'astot', 'romix', 'roc', 'rsurf', 'romix+rsurf'],
                  par='romix', sensor = None, add_interface = False, add_dutott = False, quiet = True):

        #######################
        ### start LUT importing
        if lut_base is None: lut_base = ac.settings['run']['lut_base']
        if lut_models is None: lut_models = ac.settings['run']['lut_models']
        if lut_pressures is None: lut_pressures = ac.settings['run']['lut_pressures']
        if lut_base_interface is None: lut_base_interface = ac.settings['run']['lut_base_interface']
        lut_dir = '{}/{}'.format(ac.config['lut_dir'], lut_base)

        ## make sure we read the required parameters
        lut_par = [p for p in lut_par_defaults]
        check_par = []
        if par not in lut_par: check_par.append(par) ## add par
        if ('romix+rsky_t' in lut_par) | (par == 'romix+rsky_t'):
            add_interface = True
            check_par += ['utott', 'dtott', 'astot', 'romix', 'rsky_t']
        if ('dutott' in lut_par) | (add_dutott):
            add_dutott = True
            check_par += ['dutott', 'utott', 'dtott']
        if ('romix+rsurf' in lut_par) | (par == 'romix+rsurf'):
            check_par += ['romix', 'rsurf']
        ## check if we need to add parameters to lut_par
        for k in check_par:
            if k not in lut_par: lut_par.append(k)

        ## reduce LUT loading if pressure range is known a priori
        if pressure_range is not None:
            if len(pressure_range) == 2:
                lut_pressures_float = np.asarray(lut_pressures, dtype = np.float32)
                imin = np.nanargmax(np.where((lut_pressures_float-pressure_range[0]) < 0, lut_pressures_float, np.nan))
                imax = np.nanargmin(np.where((lut_pressures_float-pressure_range[1]) > 0, lut_pressures_float, np.nan))
                lut_pressures = [lut_pressures[i] for i in range(imin, imax+1)]

        ## run through models and pressures for loading LUTs
        self.lutd = {}
        for m in lut_models:
            ## load aerlut for each pressure
            lut = '{}-MOD{}'.format(lut_base, m)
            for p in lut_pressures:
                lutid = '{}-{}mb'.format(lut, '{}'.format(p).zfill(4))
                if not quiet: print('Importing LUT: MOD{}, pressure {}, {}'.format(m, p, lutid))
                lut_data_, lut_meta_ = ac.aerlut.import_lut(lutid, lut_dir, sensor = sensor, lut_par = None, get_remote = get_remote)
                if p == lut_pressures[0]:
                    if sensor is None:
                        lut_data = [lut_data_]
                    else: ## to stack per band
                        lut_data = {b: [lut_data_[b]] for b in lut_data_}
                    lut_meta = {k: lut_meta_[k] for k in lut_meta_}
                    lut_meta['base'] = lut
                    lut_meta['press'] = [lut_meta['press']]
                else:
                    if sensor is None:
                        lut_data.append(lut_data_)
                    else: ## to stack per band
                        for b in lut_data: lut_data[b].append(lut_data_[b])
                    lut_meta['press'].append(lut_meta_['press'])
                del lut_data_
                del lut_meta_

            ## stack lut
            if sensor is None:
                lut_data = np.stack(lut_data)
            else:
                for b in lut_data: lut_data[b] = np.stack(lut_data[b])
            ## generate parameter index
            ipd = {p:i for i,p in enumerate(lut_meta['par'])}

            ## add interface if needed
            if ('romix+rsky_t' in lut_par):
                ## only at one pressure
                if not quiet: print('Importing interface LUT: MOD{}, {}'.format(m, lut_base_interface))
                int_lut_data, int_lut_meta = ac.aerlut.import_interface_lut(m, lut_base_interface = lut_base_interface,
                                                                            sensor=sensor, get_remote = get_remote)
                ## update wind meta
                lut_meta['wnd'] = int_lut_meta['wind']
                lut_meta['par'] += ['rsky_s', 'rsky_t', 'romix+rsky_t'] ## additional parameters
                ipd = {p:i for i,p in enumerate(lut_meta['par'])} ## update ipd
                ax = len(lut_meta['par'])-3 ## take into account we are adding three parameters

                ## update luts
                if sensor is None: ## generic
                    ## repeat lut for winds
                    tlut = np.repeat(lut_data, len(lut_meta['wnd']), axis=-2)
                    ## repeat for pressures
                    rlut = np.repeat(int_lut_data[np.newaxis,:], len(lut_meta['press']), axis=0)
                    del int_lut_data
                    ## insert into lut
                    tlut = np.insert(tlut, (ax), rlut, axis=1)
                    del rlut
                    ## model interface reflectance at toa
                    ## (utott * dtott * rint) / (1. - rint * astot)
                    tmp = (tlut[:, ipd['utott'],:,:,:,:,:,:]*\
                           tlut[:, ipd['dtott'],:,:,:,:,:,:]*
                           tlut[:, ax,:,:,:,:,:,:]) /\
                           (1.-tlut[:, ax,:,:,:,:,:,:] *\
                           tlut[:, ipd['astot'],:,:,:,:,:,:])
                    tlut = np.insert(tlut, (ax+1), tmp, axis=1)
                    del tmp
                    ## add romix+rsky
                    i4 = {p:i for i,p in enumerate(lut_meta['par'])}['romix']
                    tmp = tlut[:, ipd['romix'],:,:,:,:,:,:] + tlut[:, ax+1,:,:,:,:,:,:]
                    lut_data = np.insert(tlut, (ax+2), tmp, axis=1)
                    del tlut, tmp
                else: ## sensor specific
                    for b in lut_data:
                        ## repeat lut for winds
                        tlut = np.repeat(lut_data[b], len(lut_meta['wnd']), axis=-2)
                        ## repeat for pressures
                        rlut = np.repeat(int_lut_data[b][np.newaxis,:], len(lut_meta['press']), axis=0)
                        ## insert into lut
                        tlut = np.insert(tlut, (ax), rlut, axis=1)
                        del rlut
                        ## model interface reflectance at toa
                        ## (utott * dtott * rint) / (1. - rint * astot)
                        tmp = (tlut[:, ipd['utott'],:,:,:,:,:]*\
                               tlut[:, ipd['dtott'],:,:,:,:,:]*
                               tlut[:, ax,:,:,:,:,:]) /\
                               (1.-tlut[:, ax,:,:,:,:,:] *\
                               tlut[:, ipd['astot'],:,:,:,:,:])
                        tlut = np.insert(tlut, (ax+1), tmp, axis=1)
                        del tmp
                        ## add romix+rsky
                        i4 = {p:i for i,p in enumerate(lut_meta['par'])}['romix']
                        tmp = tlut[:, ipd['romix'],:,:,:,:,:] + tlut[:, ax+1,:,:,:,:,:]
                        lut_data[b] = np.insert(tlut, (ax+2), tmp, axis=1)
                        del tlut, tmp
                    del int_lut_data
            ## end add interface

            ## add romix+rsurf if needed
            if ('romix+rsurf' in lut_par):
                lut_meta['par'] += ['romix+rsurf'] ## add additional parameter
                ipd = {p:i for i,p in enumerate(lut_meta['par'])} ## update ipd
                ax = len(lut_meta['par']) -1 ## append at the end
                ## compute sum and insert in array
                if sensor is None:
                    tmp = lut_data[:,ipd['romix'],:,:,:,:,:,:] + lut_data[:,ipd['rsurf'],:,:,:,:,:,:]
                    lut_data = np.insert(lut_data, (ax), tmp, axis=1)
                    del tmp
                else:
                    for b in lut_data:
                        tmp = lut_data[b][:,ipd['romix'],:,:,:,:,:] + lut_data[b][:,ipd['rsurf'],:,:,:,:,:]
                        lut_data[b] = np.insert(lut_data[b], (ax), tmp, axis=1)
                        del tmp
            ## end add romix+rsurf

            ## add dutott if needed
            if ('dutott' in lut_par):
                lut_meta['par'] += ['dutott'] ## add additional parameter
                ipd = {p:i for i,p in enumerate(lut_meta['par'])} ## update ipd
                ax = len(lut_meta['par']) -1 ## append at the end
                ## compute product and insert in array
                if sensor is None:
                    tmp = lut_data[:,ipd['utott'],:,:,:,:,:,:] * lut_data[:,ipd['dtott'],:,:,:,:,:,:]
                    lut_data = np.insert(lut_data, (ax), tmp, axis=1)
                    del tmp
                else:
                    for b in lut_data:
                        tmp = lut_data[b][:,ipd['utott'],:,:,:,:,:] * lut_data[b][:,ipd['dtott'],:,:,:,:,:]
                        lut_data[b] = np.insert(lut_data[b], (ax), tmp, axis=1)
                        del tmp
            ## end add dutott if needed

            ## add result to lutd
            self.lutd[lut] = {'lut': lut_data, 'meta': lut_meta, 'ipd': ipd}
        ### end LUT importing
        #######################

        ## set lut info
        self.luts = list(self.lutd.keys())
        self.models = {k.split('-MOD')[1]: k for k in self.luts}

        ## get sensor information
        self.sensor = sensor
        self.bands = None
        if self.sensor is not None: ## could test whether sensor is one of the hyperspectral sensors
            self.bands = list(self.lutd[self.luts[0]]['lut'].keys())
            self.generic = False
        else:
            self.generic = True

        ## set up interpolator(s)
        self.lut = {}
        for model in self.models:
            lut = self.models[model]
            self.lut[model] = {}
            ## set model parameter indices
            self.lut[model]['ipd'] = self.lutd[lut]['ipd']

            ## set up model dimensions
            self.lut[model]['dim_names'] = ['pressure', 'parameter', 'raa', 'vza', 'sza', 'wind', 'aod']
            self.lut[model]['dim'] = [np.atleast_1d(self.lutd[lut]['meta']['press']),
                                      np.atleast_1d([i for i,p in enumerate(self.lutd[lut]['ipd'].keys())]), ## parameters
                                      np.atleast_1d(self.lutd[lut]['meta']['azi']),
                                      np.atleast_1d(self.lutd[lut]['meta']['thv']),
                                      np.atleast_1d(self.lutd[lut]['meta']['ths']),
                                      np.atleast_1d(self.lutd[lut]['meta']['wnd']),
                                      np.atleast_1d(self.lutd[lut]['meta']['tau']),]
            self.lut[model]['ndim'] = len(self.lut[model]['dim'])

            ## generic lut
            if self.generic:
                self.lut[model]['wavelengths'] = np.atleast_1d(self.lutd[lut]['meta']['wave'])
                ## set up generic model dimensions - bit of a workaround
                self.lut[model]['dim_names_generic'] = ['pressure', 'parameter', 'wavelength', 'raa', 'vza', 'sza', 'wind', 'aod']
                self.lut[model]['dim_generic'] = [np.atleast_1d(self.lutd[lut]['meta']['press']),
                                           np.atleast_1d([i for i,p in enumerate(self.lutd[lut]['ipd'].keys())]), ## parameters
                                           np.atleast_1d(self.lutd[lut]['meta']['wave']),
                                           np.atleast_1d(self.lutd[lut]['meta']['azi']),
                                           np.atleast_1d(self.lutd[lut]['meta']['thv']),
                                           np.atleast_1d(self.lutd[lut]['meta']['ths']),
                                           np.atleast_1d(self.lutd[lut]['meta']['wnd']),
                                           np.atleast_1d(self.lutd[lut]['meta']['tau']),]
                self.lut[model]['ndim_generic'] = len(self.lut[model]['dim_generic'])
                ## set up interpolator
                self.lut[model]['rgi'] = scipy.interpolate.RegularGridInterpolator(self.lut[model]['dim_generic'],
                                                                                   self.lutd[lut]['lut'][:,:,:,:,:,:,:,:],
                                                                                   bounds_error=False, fill_value=np.nan)
            ## sensor specific model
            else:
                self.lut[model]['rgi'] = {}
                for bi, b in enumerate(self.bands):
                    self.lut[model]['rgi'][b] = scipy.interpolate.RegularGridInterpolator(self.lut[model]['dim'],
                                                                                          self.lutd[lut]['lut'][b][:,:,:,:,:,:,:],
                                                                                          bounds_error=False, fill_value=np.nan)

            ## make dim dict
            self.lut[model]['dim_dict'] = {d: self.lut[model]['dim'][di] for di, d in enumerate(self.lut[model]['dim_names'])}

        ## check if we keep the lut data in memory
        self.keep_lut = keep_lut
        if not keep_lut: del self.lutd

    ## do the interpolation
    def __call__(self, model, parameter, xi, bands = None, wavelengths = None):

        ## identify lut
        ## lut is the name of the lut, mi the aerosol model index
        lut = None
        if model in self.lut:
            mi = '{}'.format(model)
            lut = self.models[mi]
        elif model in self.luts:
            for li, mi in enumerate(self.models):
                if self.models[mi] == model:
                    lut = self.models[mi]
        else:
            for li, mi in enumerate(self.models):
                if '{}'.format(model) == '{}'.format(mi):
                    lut = self.models[mi]
                    break
        if lut is None:
            print('Could not identify model {}'.format(model))
            print('Available models: {}'.format(list(self.models.keys())))
            return

        ## parameter is now provided as separate argument, but could with some modification also be passed in tuple or dict
        ## test lut dimensions (add one for parameter which is provided separately as string here)
        nxi = len(xi)
        if parameter is not None: nxi += 1
        #if (wavelengths is None) & (self.generic): nxi +=1
        if nxi != self.lut[mi]['ndim']:
            print('Number of dimensions ({}) does not match with those in LUT ({})'.format(nxi,  self.lut[mi]['ndim']))
            return

        ## identify index of requested parameter (passed as string)
        par = None
        if parameter in self.lut[mi]['ipd']:
            par = self.lut[mi]['ipd'][parameter]
        else:
            print('Could not find parameter {} in lut {} (model {})'.format(parameter, lut, mi))
            return

        ## test dimensions of requested interpolation
        if type(xi) in (list, tuple): ## assume dimensions are in order
            isize = [np.atleast_1d(x).size for x in xi]
            ishape = [np.atleast_1d(x).shape for x in xi]
        elif type(xi) == dict: ## set dimensions based on dim_names
            isize = []
            ishape = []
            for di, dim in enumerate(self.lut[model]['dim_names']):
                if dim == 'parameter':
                    continue
                elif dim == 'wavelength':
                    continue
                else:
                    isize += [np.atleast_1d(xi[dim]).size]
                    ishape += [np.atleast_1d(xi[dim]).shape]
        else:
            print('Requested interpolation needs to be list, tuple or dict')
            return

        ## check if we can understand the requested interpolation
        ssize = set(isize)
        if ((len(ssize)!=1) & (1 not in ssize)) | (len(ssize)>2):
            print('Inconsistent dimensions for requested interpolation: {}'.format(isize))
            return
        if ((len(ssize)!=1) & (self.generic)):
            print('Currently only single point interpolations supported for generic LUT')
            print('Dimensions for requested interpolation: {}'.format(isize))
            return

        ## set up wavelength dimension
        ## in generic case the dim and dim_names exclude wavelengths, but the dimension is there in the rgi (and in dim_generic)
        if self.generic:
            if wavelengths is None:
                wave = np.asarray(self.lut[model]['wavelengths']).flatten()
            else:
                wave = np.asarray(wavelenths).flatten()
            wave = wave.T

        ## find max size and get shape
        xil = max(isize)
        ishape = [ishape[i] for i, x in enumerate(isize) if x == xil][0]
        ## add parameter index to interpolator
        pard = [par]
        if xil > 1: pard = [[par] * xil]
        if type(xi) in (list, tuple): ## assume dimensions are in order
            xi_ = [xi[0]] + pard + list(xi[1:])
            #if not self.generic:
            #    xi_ = [xi[0]] + pard + list(xi[1:])
            #else:
            #    xi_ = [xi[0]] + pard + wave + list(xi[1:])
        elif type(xi) == dict: ## set dimensions based on dim_names
            xi_ = []
            for di, dim in enumerate(self.lut[model]['dim_names']):
                if dim == 'parameter':
                    xi_ += pard
                elif dim == 'wavelength':
                    continue
                    #if not self.generic: continue
                    #xi_ += wave
                else:
                    xi_ += [xi[dim]]

        ## unify dimensions
        for di, dim in enumerate(self.lut[model]['dim_names']):
            ## make array so it can be filled
            xi_[di] = np.atleast_1d(xi_[di])#.flatten()
            ## if any LUT dimension == 1, set the requested position at the only LUT point (otherwise rgi will give nans)
            if len(self.lut[model]['dim'][di]) == 1: xi_[di][:] = self.lut[model]['dim'][di][0]
            ## repeat elements if the current requested length is not the longest
            if xi_[di].size != xil: xi_[di] = np.repeat(xi_[di], xil)
            ## reshape to the first size
            xi_[di] = xi_[di].reshape(ishape)

        ## convert to array and transpose to get correct shape of interpolation request
        xi_ = np.asarray(xi_).T

        ## run through bands and do interpolation
        if (self.generic) | (self.bands is None):
            #print((xi_[0,0], xi_[0,1], wave, xi_[0,2], xi_[0,3], xi_[0,4], xi_[0,5],xi_[0,6]))
            ret = self.lut[mi]['rgi']((xi_[0,0], xi_[0,1], wave, xi_[0,2], xi_[0,3], xi_[0,4], xi_[0,5],xi_[0,6]))
            return(ret)
        else:
            ## list interpolation bands - if bands is None return all bands
            ibands = [b for b in self.bands]
            if bands is not None:
                if type(bands) is not list: bands = [bands]
                for b in bands:
                    if b not in self.bands:
                        print('Band {} not in sensor specific LUT: {}'.format(b, self.sensor))
                        print('Available bands: {}'.format(self.bands))
                        return
                ibands = [b for b in bands]

            ## interpolate for bands
            ret = [self.lut[mi]['rgi'][b]((xi_))[0:xil] for b in ibands]
            if len(ibands) == 1: ret = ret[0] ## 0 index if only one band is requested
            ret = np.asarray(ret)
            return(ret)
