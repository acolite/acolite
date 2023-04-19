## QV Jul 2019
## import libradtran simulation results
def read_sim(sfil):
    sdata = None

    with open(sfil, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line[0] in ['#']: continue
            sp = line.split(',')
            if sdata is None:
                sheader = sp
                sdata = {h:[] for h in sheader}
            else:
                for ih, h in enumerate(sheader):
                    v = sp[ih]
                    if ih > 0: v = float(v)
                    sdata[h].append(v)
    return(sdata)
