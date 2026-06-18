## hydrolight.dfile.read
## reads Hydrolight outputfile
##
## written by Quinten Vanhellemont, RBINS
## 2018-04-23
## modifications: 2022-11-27 (QV) added ECOLIGHT
##                2026-06-16 (QV) new function, renamed wavelengths to wavelength

def read(file):

    data = {'wavelength':[]}
    wave_data = {}
    with open(file, 'r') as f:
        for li, line in enumerate(f.readlines()):
            line = line.strip()

            if line[0:10] == 'HYDROLIGHT':
                if 'title' not in data:
                    data['software'] = 'HYDROLIGHT'
                    data['title'] = line
                continue

            if line[0:8] == 'ECOLIGHT':
                if 'title' not in data:
                    data['software'] = 'ECOLIGHT'
                    data['title'] = line
                continue

            if line[0:4] == 'NOTE':
                if 'note' not in data:
                    data['note'] = line
                continue

            if line[0:3] == 'nmu': #li == 2:
                split = line.split('=')
                grid_h = split[0].split(',')
                grid_v = split[1].split()
                for gi, g in enumerate(grid_h):
                    gn = g.strip()
                    if gn not in data: data[gn] = int(grid_v[gi])
                continue

            if line[0:10] == 'wavelength':
                split = line.split()
                try:
                    cur_band = int(split[2])
                except:
                    cur_band = int(split[1].replace('band', ''))
                cur_wave = float(split[-1])
                if cur_wave not in data:
                    data[cur_wave] = {}
                    data['wavelength'].append(cur_wave)
                continue

            if line[0].isalpha():
                cur_par = line.strip()
                cur_parname = cur_par.split()[0]
                if cur_parname not in data[cur_wave]:
                    data[cur_wave][cur_parname] = {'title':cur_par, 'par': cur_parname, 'data': []}
                continue

            cur_data = line.split()
            if cur_par == 'imisc':
                cur_data = [int(v) for v in cur_data]
            else:
                cur_data = [float(v) for v in cur_data]

            data[cur_wave][cur_parname]['data'] += cur_data

    return(data)
