## class lira
## class to run libRadtran and read outputs
## written by Quinten Vanhellemont, RBINS
## 2024-10-30
## modifications:

import os, subprocess, time, datetime, string
import acolite as ac
import numpy as np

class lira(object):
    def __init__(self, run_name = None, output = None, libradtran_dir = None,
                override = True, quiet = True, delete = True):
        self.override = override ## overwrite inp and out files if they exist
        self.quiet = quiet ## do not printout information
        self.delete = delete # delete inp and out files

        self.working_dir = os.getcwd()

        ## get libradtran dir - used for data path and running bin/uvspec
        if libradtran_dir is not None:
            self.libradtran_dir = libradtran_dir
        else:
            self.libradtran_dir = ac.config['libradtran_dir']

        ## get output dir - use acolite/scratch dir if not set
        if output is not None:
            self.output = output
        else:
            self.output = ac.config['scratch_dir'] + os.sep + 'libRadtran'

        ## get run name - use datetime if not set
        if run_name is not None:
            self.run_name = run_name
        else:
            self.run_name = datetime.datetime.now().strftime('%Y%m%d_%H%M%S') + \
                            '_{}'.format(''.join(np.random.choice(list(string.ascii_uppercase), size=6, replace=True)))

        self.set_file_paths()

        self.cfg = {}
        self.set_defaults()

        self.result = None

    def __del__(self):
        if self.delete:
            if not self.quiet: print('Deleting simulation {}'.format(self.run_name))
            self.clear_inp_file()
            self.clear_out_file()

    def set_file_paths(self):
        self.inp_file = self.output + os.sep + self.run_name + '.INP.txt'
        self.out_file = self.output + os.sep + self.run_name + '.OUT.txt'

    def clear_inp_file(self):
        if os.path.exists(self.inp_file):
            if not self.quiet: print('Deleting {}'.format(self.inp_file))
            os.remove(self.inp_file)

    def clear_out_file(self):
        if os.path.exists(self.out_file):
            if not self.quiet: print('Deleting {}'.format(self.out_file))
            os.remove(self.out_file)

    def set_defaults(self):
        self.cfg["data_files_path"] =  '{}/data'.format(self.libradtran_dir)
        self.cfg["output_user"] = 'lambda eglo edir edn eup uu'
        self.cfg["quiet"] = ''
        self.cfg["mol_modify O3"] = '300 DU'
        self.cfg["mol_modify H2O"] = '15 MM'
        self.cfg["albedo"] = '0'
        self.cfg["zout" ] = 'SUR TOA'
        self.cfg["pressure"] = '1013.25'

        self.cfg["mol_abs_param reptran"] = 'coarse'
        self.cfg["rte_solver"] = 'disort'
        self.cfg["number_of_streams"] = '16'

        self.cfg["sza"] = '30'
        self.cfg["phi0"] = '0'
        self.cfg["wavelength"] = '300 2500'
        self.cfg["umu"] = '1.0' ## cos view angles
        self.cfg["phi"] = '0' ## view azimuths

    def write_inp_file(self):
        self.set_file_paths() ## set file paths in case they changed
        if not os.path.exists(os.path.dirname(self.inp_file)):  os.makedirs(os.path.dirname(self.inp_file))
        ## write key setting combinations
        with open(self.inp_file, 'w', encoding='utf-8') as f:
            for s in self.cfg: f.write('{} {}\n'.format(s, self.cfg[s]))

    def read_out_file(self):
        data = np.loadtxt(self.out_file)
        self.result = {}
        levels = self.cfg['zout'].split()
        pars = self.cfg['output_user'].split()
        for i, l in enumerate(levels):
            self.result[l] = {}
            for j, p in enumerate(pars):
                self.result[l][p] = data[i::2,j]
            ## add total irradiance if not present
            if ('eglo' not in self.result[l]) & \
               ('edir' in self.result[l]) & ('edn' in self.result[l]):
               self.result[l]['eglo'] = self.result[l]['edir'] + self.result[l]['edn']
        del data

    def run(self):
        ## change parameters for twostr solver if still defaults
        if self.cfg["rte_solver"] == 'twostr':
            if (self.cfg["output_user"] == 'lambda eglo edir edn eup uu'):
                self.cfg["output_user"] = 'lambda eglo edir edn eup uavg'

        if self.override:
            self.clear_inp_file()
            self.clear_out_file()

        ## write settings
        if not os.path.exists(self.inp_file): self.write_inp_file()

        ## clear empty output file
        if os.path.exists(self.out_file):
            if (os.stat(self.out_file).st_size == 0): self.clear_out_file()

        ## run simulation
        if not os.path.exists(self.out_file):
            t0 = time.time()
            if not self.quiet: print('Running simulation {}'.format(self.run_name))
            ## change directory to libradtran and run uvspec
            os.chdir(self.libradtran_dir+'/bin')
            cmd = ['./uvspec','< {}'.format(self.inp_file),'> {}'.format(self.out_file)]
            p = subprocess.run(' '.join(cmd), shell=True, stdout=subprocess.PIPE)
            os.chdir(self.working_dir)
            t1 = time.time()
            if not self.quiet:
                print('Finished simulation {}'.format(self.run_name))
                print('Simulation took {:.1f}s'.format(t1-t0))

        ## read output
        if os.path.exists(self.out_file):
            if (os.stat(self.out_file).st_size > 0):
                if not self.quiet: print('Reading simulation {}'.format(self.run_name))
                self.read_out_file()
