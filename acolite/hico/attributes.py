## get some attributes from HICO L1B NetCDF file
## QV 2021-08-03

def attributes(file):
    import h5py
    atts = {}
    with h5py.File(file, mode='r') as f:
        df = f['metadata']['FGDC']['Identification_Information']['Time_Period_of_Content']
        for k in df.attrs.keys(): atts[k] = df.attrs[k].decode('utf-8')
        df = f['metadata']['FGDC']['Identification_Information']['Platform_and_Instrument_Identification']
        for k in df.attrs.keys(): atts[k] = df.attrs[k].decode('utf-8')
    return(atts)
