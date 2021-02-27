## QV Jul 2019

def read_out(outfile, parameters = ['lambda','eup','uu']):
    import numpy as np
    
    data = {}
    for d in ['SUR', 'TOA']: 
        data[d] = {par:[] for par in parameters}
        
    with open(outfile) as f:
        for il, line in enumerate(f.readlines()):
            line = line.strip()
            sp = [float(v) for v in line.split()]

            if il % 2 == 0: tag = 'SUR'
            if il % 2 == 1: tag = 'TOA'
            
            for ip, par in enumerate(parameters):
                data[tag][par].append(sp[ip])
                
    for k in data:
        for l in data[k]:
            data[k][l] = np.asarray(data[k][l])
    return(data)
