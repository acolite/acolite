### hydrolight.write
## write simple HE5.2/5.3 parameter file
##
## written by Quinten Vanhellemont, RBINS
## 2022-12-02
## modifications: 2026-06-16 (QV) new function

def write(file, data, comment='custom HL file'):
    import os
    header = ['blank']*9
    header[0] = comment
    name = [k for k in data.keys() if k != 'waves'][0]

    if not os.path.exists(os.path.dirname(file)):
        os.makedirs(os.path.dirname(file))

    with open(file, 'w', encoding='utf-8') as f:
        ## write blank header
        for i, c in enumerate(header):
            f.write(c)
            f.write('\n')
        ## write column names
        f.write('\t'.join(['wavelen(nm)', name]))
        f.write('\n')

        ## write data
        for i, w in enumerate(data['waves']):
            if w == -1: continue
            f.write('\t'.join([str(s) for s in [w, data[name][i]]]))
            f.write('\n')
        f.write('\t'.join([str(s) for s in [-1, -1.0, -1.0]]))
        f.write('\n')
