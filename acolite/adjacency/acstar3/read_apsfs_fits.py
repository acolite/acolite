## reads apsfs cumulative integral PSF fits
## QV 2021-07-07
## modifications:   2021-10-06 (QV) updates for new fit files

def read_apsfs_fits(model_name, sensor):
    import acolite as ac
    import os

    #rayf = ac.config['data_dir']+'/Shared/apsfs_model/fa_rayleigh.csv'
    #aerf = ac.config['data_dir']+'/Shared/apsfs_model/fa_{}_{}.csv'.format(sensor.split('_')[1].lower(), model_name[0:3])
    ss = sensor.lower().split('_')
    if sensor in ['L8_OLI']:
        rayf = ac.config['data_dir']+'/Shared/psf_fit/{}_{}_R_res030_v00_p1100_annular_coeff.csv'.format(ss[1], ss[0])
        aerf = ac.config['data_dir']+'/Shared/psf_fit/{}_{}_{}_res030_v00_p1013_annular_coeff.csv'.format(ss[1], ss[0], model_name[0].upper())
        band_names = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    if sensor in ['S2A_MSI', 'S2B_MSI']:
        rayf = ac.config['data_dir']+'/Shared/psf_fit/{}_{}_R_res010_v00_p1100_annular_coeff.csv'.format(ss[1], ss[0])
        aerf = ac.config['data_dir']+'/Shared/psf_fit/{}_{}_{}_res010_v00_p1013_annular_coeff.csv'.format(ss[1], ss[0], model_name[0].upper())
        band_names = ['1', '2', '3', '4', '5', '6', '7', '8', '8A', '9', '10', '11', '12']
    if 'PlanetScope' in sensor:
        sfam = sensor[-2:]
        if sfam in ['0c','0d']:
            sfam = '0c0d'
        if sfam in ['0f','10']:
            sfam = '0f10'
        rayf = ac.config['data_dir']+'/Shared/psf_fit/{}_{}_R_res003_v00_p1100_annular_coeff.csv'.format(sfam, 'ps')
        aerf = ac.config['data_dir']+'/Shared/psf_fit/{}_{}_{}_res003_v00_p1013_annular_coeff.csv'.format(sfam, 'ps', model_name[0].upper())
        band_names = ['Blue', 'Green', 'Red', 'NIR']

    print(rayf, aerf)

    rdata, adata = {}, {}
    if os.path.exists(rayf):
        with open(rayf, 'r', encoding='utf-8') as f:
            for il, line in enumerate(f.readlines()):
                line = line.strip()
                sp = [s.strip().strip('"') for s in line.split(',')]
                if il == 0:
                    header = [s.replace(' ', '') for s in sp[1:]]
                else:
                    rdata = {h: float(sp[ih+1]) for ih, h in enumerate(header)}
    if os.path.exists(aerf):
        with open(aerf, 'r', encoding='utf-8') as f:
            for il, line in enumerate(f.readlines()):
                line = line.strip()
                sp = [s.strip().strip('"') for s in line.split(',')]
                if il == 0:
                    header = [s.replace(' ', '') for s in sp[1:]]
                else:
                    #adata[sp[0].strip('band_')] = {h: float(sp[ih+1]) for ih, h in enumerate(header)}
                    adata[band_names[int(sp[0].strip('band_'))-1]] = {h: float(sp[ih+1]) for ih, h in enumerate(header)}
    return(rdata, adata)
