## QV 2020-07-24
## QV 2021-12-22 adapted for MERIS

def meris_band_info():
    import acolite as ac

    band_file = '{}/{}/{}'.format(ac.config['data_dir'], 'EN1', 'band_info_meris.txt')
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

                band_data['M{}'.format(str(c['band']).zfill(2))] = c
    return(band_data)
