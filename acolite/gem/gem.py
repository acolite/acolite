## simple NetCDF gem object
## written by Quinten Vanhellemont, RBINS
## 2021-04-01
## modifications: 2021-04-01 (QV) added some write support
##                2021-12-08 (QV) added nc_projection
##                2022-02-15 (QV) added L9/TIRS

import acolite as ac
import os, sys
import numpy as np
from netCDF4 import Dataset

class gem(object):
        def __init__(self, file, new=False,
                    netcdf_compression=False,
                    netcdf_compression_level=4,
                    netcdf_compression_least_significant_digit=None):
            self.file=file
            self.data_mem = {}
            self.data_att = {}
            self.store = False
            self.bands = {}
            self.verbosity = 0
            self.nc_projection = None

            self.netcdf_compression=netcdf_compression
            self.netcdf_compression_level=netcdf_compression_level
            self.netcdf_compression_least_significant_digit=netcdf_compression_least_significant_digit

            if new:
                self.new = True
                if os.path.exists(self.file):
                    os.remove(self.file)

            if os.path.exists(self.file):
                self.new = False
                self.gatts_read()
                self.datasets_read()
                self.nc_projection = ac.shared.nc_read_projection(self.file)

        def gatts_read(self):
            self.gatts = ac.shared.nc_gatts(self.file)
            ## detect thermal sensor
            if self.gatts['sensor'] == 'L8_OLI':
                self.gatts['thermal_sensor'] = 'L8_TIRS'
                self.gatts['thermal_bands'] = ['10', '11']
            elif self.gatts['sensor'] == 'L9_OLI':
                self.gatts['thermal_sensor'] = 'L9_TIRS'
                self.gatts['thermal_bands'] = ['10', '11']
            elif self.gatts['sensor'] == 'L5_TM':
                self.gatts['thermal_sensor'] = 'L5_TM'
                self.gatts['thermal_bands'] = ['6']
            elif self.gatts['sensor'] == 'L7_ETM':
                self.gatts['thermal_sensor'] = 'L7_ETM'
                self.gatts['thermal_bands'] = ['6_vcid_1', '6_vcid_2']

        def datasets_read(self):
            self.datasets = ac.shared.nc_datasets(self.file)

        def data(self, ds, attributes=False, store=False, return_data=True):
            if ds in self.data_mem:
                cdata = self.data_mem[ds]
                if ds in self.data_att:
                    catt = self.data_att[ds]
                else:
                    catt = {}
            else:
                if ds in self.datasets:
                    cdata, catt = ac.shared.nc_data(self.file, ds, attributes=True)
                    cmask = cdata.mask
                    cdata = cdata.data
                    if cdata.dtype in [np.dtype('float32'), np.dtype('float64')]:
                        cdata[cmask] = np.nan
                    if (self.store) or (store):
                        self.data_mem[ds] = cdata
                        self.data_att[ds] = catt
                else:
                    return()
            if return_data:
                if attributes:
                    return(cdata, catt)
                else:
                    return(cdata)

        def write(self, ds, data, ds_att = {}):
            if self.new:
                if os.path.exists(self.file):
                    os.remove(self.file)
            ac.output.nc_write(self.file, ds, data, attributes=self.gatts,
                                dataset_attributes=ds_att, new=self.new,
                                nc_projection=self.nc_projection,
                                netcdf_compression=self.netcdf_compression,
                                netcdf_compression_level=self.netcdf_compression_level,
                                netcdf_compression_least_significant_digit=self.netcdf_compression_least_significant_digit)
            if self.verbosity > 0: print('Wrote {}'.format(ds))
            self.new = False

        def update_attributes(self):
            with Dataset(self.file, 'a', format='NETCDF4') as nc:
                for key in self.gatts.keys():
                    if self.gatts[key] is not None:
                        try:
                            setattr(nc, key, self.gatts[key])
                        except:
                            print('Failed to write attribute: {}'.format(key))
