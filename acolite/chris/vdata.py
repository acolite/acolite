## read CHRIS vdata Gains and Mode data
##
##
## QV 2018-03

def vdata(file):
    from pyhdf import HDF
    from pyhdf import VS
    f = HDF.HDF(file)
    vs = f.vstart()                  # init vdata interface

    ## read Gains table
    vd = vs.attach('Gain Information')
    fields = vd._fields
    gains = {}

    while 1:
        try:
            rec = vd.read()[0]
            gains[rec[0]] = float(rec[1])
        except:
            break
    vd.detach()

    ## read Modes table
    vd = vs.attach('Mode Information')      # attach 'INVENTORY' in read mode
    fields = vd._fields

    mode_info = {}
    r=0
    while 1:
        try:
            rec = vd.read()[0]
            cd = {fn:rec[fi] for fi,fn in enumerate(fields)}
            mode_info[r] = cd
            r+=1
        except:
            break
    vd.detach()

    ## close vdata
    vs.end()                  # terminate the vdata interface
    f.close()                 # close the HDF file

    return(gains, mode_info)
