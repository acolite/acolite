## QV 2020-07-24

def olci_band_info():
    import acolite as ac

    band_file = '{}/{}/{}'.format(ac.config['data_dir'],'S3', 'band_info_olci.txt')
    band_data = {}
    with open(band_file, 'r') as f:
        for il, line in enumerate(f.readlines()):
            line = line.strip()
            sp = line.split()
            if il == 0:
                header = sp
            else:
                c = {h:sp[ih] for ih, h in enumerate(header)}
                c['wavelength'] = float(c['wavelength']) ## changed from OBPG "reflectance" naming
                c['E0'] = float(c['E0'])
                for t in c:
                    if type(c[t]) == str: c[t] = int(c[t])

                band_data['Oa{}'.format(str(c['band']).zfill(2))] = c
    return(band_data)
