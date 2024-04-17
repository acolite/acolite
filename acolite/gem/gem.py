## simple NetCDF gem object
## written by Quinten Vanhellemont, RBINS
## 2021-04-01
## modifications: 2021-04-01 (QV) added some write support
##                2021-12-08 (QV) added nc_projection
##                2022-02-15 (QV) added L9/TIRS
##                2023-07-12 (QV) removed netcdf_compression settings
##                2024-01-31 (QV) added skip_attributes
##                2024-04-16 (QV) moved Landsat thermals to file
##                                gem is now independent of shared.nc_read
##                                file handle is passed to output.nc_write if file exists

import acolite as ac
import os, sys, json
import numpy as np
from netCDF4 import Dataset

class gem(object):
        def __init__(self, file, new = False, verbosity = 0):
            self.file=file
            self.data_mem = {}
            self.data_att = {}
            self.store = False
            self.bands = {}
            self.verbosity = verbosity

            ## nc handle and mode
            self.nc = None
            self.nc_mode = None
            self.new = new

            ## attributes and datasets
            self.gatts = None
            self.datasets = []
            self.nc_projection = None
            self.nc_projection_keys = []

            ## if creating new dataset, delete existing file
            if (self.new) & (os.path.exists(self.file)): os.remove(self.file)

            ## if exists, read in gatts, datasets, and projection
            if os.path.exists(self.file): self.setup()

        ## basic set up read in gatts, datasets, and projection
        def setup(self):
            self.new = False
            self.gatts_read()
            self.datasets_read()
            self.projection_read()

        ## create nc file handle
        def open(self, *h):
            ## set read/write mode, defaults to r
            mode = 'r' if len(h) == 0 else h[0]

            ## close file if open in wrong mode
            if (self.nc is not None) & (self.nc_mode != mode):
                if self.verbosity > 5: print('File {} already open in mode {}'.format(self.file, self.nc_mode))
                self.close()

            ## open in requested mode
            if (self.nc is None):
                if self.verbosity > 5: print('Opening {} in mode {}'.format(self.file, mode))
                self.nc = Dataset(self.file, mode)
                self.nc_mode = mode

        ## close nc
        def close(self):
            if (self.nc is not None):
                self.nc.close()
                self.nc = None
                self.nc_mode = None
                if self.verbosity > 5: print('Closed {}'.format(self.file))

        ## read attributes
        def gatts_read(self):
            if self.nc_mode != 'r': self.open('r')
            if self.gatts is None:
                self.gatts = {attr : getattr(self.nc,attr) for attr in self.nc.ncattrs()}
                ## find out if we have Landsat thermal bands
                if 'sensor' in self.gatts:
                    with open(ac.config['data_dir'] + '/Landsat/thermal_sensor.json', 'r', encoding = 'utf-8') as f:
                        thermal_dict = json.loads(f.read())
                    if self.gatts['sensor'] in thermal_dict:
                        for k in thermal_dict[self.gatts['sensor']]:
                            self.gatts[k] = thermal_dict[self.gatts['sensor']][k]

        ## read available datasets
        def datasets_read(self, group = None):
            if self.nc_mode != 'r': self.open('r')
            if group is not None:
                if group in self.nc.groups: self.datasets = list(self.nc.groups[group].variables.keys())
            else:
                self.datasets = list(self.nc.variables.keys())

        ## read dataset
        def data(self, ds, attributes = False, store = False, return_data = True, sub = None):
            ## data already in memory
            if ds in self.data_mem:
                cdata = self.data_mem[ds]
                catt = {}
                if ds in self.data_att: catt = self.data_att[ds]
            ## read in dataset
            else:
                if self.nc_mode != 'r': self.open('r')
                if ds in self.datasets:
                    ## get data
                    if (ds not in self.nc_projection_keys) & (sub is not None):
                        cdata = self.nc.variables[ds][sub[1]:sub[1]+sub[3]:1,sub[0]:sub[0]+sub[2]:1]
                    else:
                        cdata = self.nc.variables[ds][:]
                    ## get attributes
                    catt = {attr : getattr(self.nc.variables[ds],attr) for attr in self.nc.variables[ds].ncattrs()}
                    ## mask data
                    cmask = cdata.mask
                    cdata = cdata.data
                    if ds not in self.nc_projection_keys:
                        if cdata.dtype in [np.dtype('float32'), np.dtype('float64')]:
                            cdata[cmask] = np.nan
                    if (self.store) or (store):
                        self.data_mem[ds] = cdata
                        self.data_att[ds] = catt
                else:
                    return

            ## return_data is set by default
            ## but this method can be used to load all datasets without returning data
            if return_data:
                if attributes:
                    return(cdata, catt)
                else:
                    return(cdata)

        ## read data for nc_projection
        def projection_read(self):
            if 'projection_key' in self.gatts:
                self.nc_projection = {}
                self.nc_projection_keys = [self.gatts['projection_key'], 'x', 'y']
                for k in self.nc_projection_keys:
                    d, a = self.data(k, attributes = True)
                    self.nc_projection[k] = {'data': d, 'attributes': a}

        ## write dataset
        def write(self, ds, data, ds_att = {}, replace_nan = False, update_projection = False):
            if self.new:
                self.close() # close if open
                if os.path.exists(self.file): os.remove(self.file) # delete if exists
                ## create new netcdf file
                self.nc = ac.output.nc_write(self.file, ds, data, dataset_attributes=ds_att,
                                             new=self.new, return_nc = True, attributes=self.gatts,
                                             nc_projection=self.nc_projection)
                self.nc_mode = 'w'
                self.new = False # new file has been created
                self.setup() # read in attributes
            else:
                if self.nc_mode != 'a': self.open('a')
                ac.output.nc_write(self.nc, ds, data, dataset_attributes = ds_att, replace_nan = replace_nan,
                                            update_projection = update_projection, nc_projection=self.nc_projection)

            if self.verbosity > 0: print('Wrote {}'.format(ds))

        ## update global attributes
        def gatts_update(self, close = True):
            if self.nc_mode != 'a': self.open('a')
            for key in self.gatts.keys():
                if key in ac.config['skip_attributes']: continue
                if self.gatts[key] is not None:
                    try:
                        setattr(self.nc, key, self.gatts[key])
                    except:
                        if self.verbosity > 3:
                            print('Failed to write attribute: {}'.format(key))
                            print('Attribute {} has type {}'.format(key, type(self.gatts[key])))
            if close: self.close() ## close to update attributes
