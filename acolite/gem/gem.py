## simple NetCDF gem object
## written by Quinten Vanhellemont, RBINS
## 2021-04-01
## modifications: 2021-04-01 (QV) added some write support
##                2021-12-08 (QV) added nc_projection
##                2022-02-15 (QV) added L9/TIRS
##                2023-07-12 (QV) removed netcdf_compression settings
##                2024-01-31 (QV) added skip_attributes
##                2024-04-16 (QV) moved Landsat thermals to file

import acolite as ac
import os, sys, json
import numpy as np
from netCDF4 import Dataset

class gem(object):
        def __init__(self, file, new=False):
            self.file=file
            self.data_mem = {}
            self.data_att = {}
            self.store = False
            self.bands = {}
            self.verbosity = 0
            self.nc_projection = None

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

            ## find out if we have Landsat thermal bands
            with open(ac.config['data_dir'] + '/Landsat/thermal_sensor.json', 'r', encoding = 'utf-8') as f:
                thermal_dict = json.loads(f.read())
            if self.gatts['sensor'] in thermal_dict:
                for k in thermal_dict[self.gatts['sensor']]:
                    self.gatts[k] = thermal_dict[self.gatts['sensor']][k]

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
                    return
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
                                nc_projection=self.nc_projection)
            if self.verbosity > 0: print('Wrote {}'.format(ds))
            self.new = False

        def update_attributes(self):
            with Dataset(self.file, 'a', format='NETCDF4') as nc:
                for key in self.gatts.keys():
                    if key in ac.config['skip_attributes']: continue
                    if self.gatts[key] is not None:
                        try:
                            setattr(nc, key, self.gatts[key])
                        except:
                            print('Failed to write attribute: {}'.format(key))
