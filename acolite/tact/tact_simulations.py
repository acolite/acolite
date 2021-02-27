## QV Jul 2019
##         2019-12-17 renamed, integrated in tact
##         2021-02-27 (QV) integrated in acolite renamed from run_thermal_sim

def tact_simulations(sonde, atmosphere="../data/atmmod/afglss.dat", obase=None,
                    pdate=None, rsr_data=None, brightness = False, override = False):
    import os
    import datetime
    import numpy as np

    import acolite as ac

    sim = os.path.basename(sonde)
    sim = os.path.splitext(sim)[0]
    sim = sim.replace('_reformatted', '')

    if pdate is None: pdate = datetime.datetime.now().strftime('%Y%m%d')
    if obase is None: obase = os.path.dirname(sonde)

    odir = '{}/{}/{}/'.format(obase, pdate, sim)
    if not os.path.exists(odir): os.makedirs(odir)

    ofile = '{}/{}_output.txt'.format(odir, sim)

    parameters = ['lambda','uu']
    data = {}
    for i in range(3):
        run = 'run{}'.format(i+1)
        if run == 'run1':
            t = -30
            e = 1.0
            look_down=True
        elif run == 'run2':
            t = 0
            e = 1.0
            look_down=True
        elif run == 'run3':
            t = -273.15
            e = 0.9
            look_down=True
        else:
            print('Run not configured')
            continue

        runfile = '{}/{}.inp'.format(odir, run)
        outfile = runfile.replace('.inp', '.out')
        if (os.path.exists(outfile)):
            statinfo = os.stat(outfile)
            if statinfo.st_size == 0:
                os.remove(outfile)

        if (override) or (os.path.exists(outfile) is False):
            sur_temperature = None
            if t is not None: sur_temperature=t+273.15

            cfg = ac.tact.libradtran_cfg(runfile=runfile,
                             look_down=look_down,
                             sur_temperature=sur_temperature,
                             atmosphere=atmosphere,
                             brightness=brightness,
                             parameters=parameters,
                             radiosonde=sonde,
                             albedo=1-e, thv=0, phi=0, phi0=0)


            outfile = ac.tact.libradtran_run(runfile)
        data[run] = ac.tact.read_out(outfile,parameters = parameters)

    ## compute Lu, tau
    waves = np.asarray(data['run1']['SUR']['lambda'])
    n = len(waves)
    muwave = waves/1000.

    #if not os.path.exists(ofile) or (override):
    if True:
        tau = []
        Lu = []

        for i in range(n):
            mxc = (data['run2']['TOA']['uu'][i] - data['run1']['TOA']['uu'][i])/\
                  (data['run2']['SUR']['uu'][i] - data['run1']['SUR']['uu'][i])


            bxc =  data['run2']['TOA']['uu'][i] - (data['run2']['SUR']['uu'][i] * mxc)
            if mxc < 0.01: bxc = np.nan

            tau.append(mxc)
            Lu.append(bxc)
        tau = np.asarray(tau)
        Lu = np.asarray(Lu)*1000

        ## compute Ld
        Ld = ((((data['run3']['TOA']['uu']*1000.)-Lu)/tau))/(1-0.9)

        Ld[np.isnan(Ld)] = 0
        Lu[np.isnan(Lu)] = 0

        ## write output
        with open(ofile, 'w') as f:
            f.write('{}\n'.format('# LibRadtran results {}'.format(sim)))
            f.write('{}\n'.format('# {}'.format(datetime.datetime.now().isoformat())))
            f.write('{}\n'.format('# contact: Quinten Vanhellemont, RBINS'))
            f.write('{}\n'.format(','.join(['wavelength,tau,Lu,Ld'])))
            for i in range(n):
                f.write('{}\n'.format(','.join([str(s) for s in [waves[i], tau[i], Lu[i], Ld[i]]])))

    if True:
        ## resample

        if rsr_data is not None:
            for satsen in rsr_data:
                ofile = '{}/{}_{}.txt'.format(odir, sim, satsen)
                if os.path.exists(ofile) & (override is False): continue
                ret_wave = ac.shared.rsr_convolute_dict(muwave,
                                                    muwave,
                                                    rsr_data[satsen]['rsr'],
                                                    wave_range=[9,14], wave_step=0.05)

                ret_tau = ac.shared.rsr_convolute_dict(muwave,
                                                     tau,
                                                     rsr_data[satsen]['rsr'],
                                                     wave_range=[9,14], wave_step=0.05)

                ret_Lu = ac.shared.rsr_convolute_dict(muwave,
                                                     Lu,
                                                     rsr_data[satsen]['rsr'],
                                                     wave_range=[9,14], wave_step=0.05)

                ret_Ld = ac.shared.rsr_convolute_dict(muwave,
                                                      Ld,
                                                      rsr_data[satsen]['rsr'],
                                                     wave_range=[9,14], wave_step=0.05)

                with open(ofile, 'w') as f:
                    f.write('{}\n'.format('# LibRadtran results {} - {}'.format(sim, satsen)))
                    f.write('{}\n'.format('# {}'.format(datetime.datetime.now().isoformat())))
                    f.write('{}\n'.format('# contact: Quinten Vanhellemont, RBINS'))
                    f.write('{}\n'.format(','.join(['band,wavelength,tau,Lu,Ld'])))
                    for b in ret_Ld:
                        f.write('{}\n'.format(','.join([str(s) for s in [b, ret_wave[b],
                                                                         ret_tau[b],
                                                                         ret_Lu[b], ret_Ld[b]]])))
