## def segment read
## reads Himawari segment file
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-19
## modifications:

def segment_read(segment_file, byte_range = None):
    import bz2
    if byte_range is not None:
        if len(byte_range) != 2:
            print('Keyword byte_range needs two elements, offset and length')
            return

    ## read binary data
    if segment_file.endswith('.bz2'):
        with bz2.open(segment_file, 'rb') as f:
            if byte_range is None:
                bin = f.read()
            else:
                f.seek(byte_range[0])
                if byte_range[1] is None:
                    bin = f.read()
                else:
                    bin = f.read(byte_range[1])

    else:
        with open(segment_file, 'rb') as f:
            if byte_range is None:
                bin = f.read()
            else:
                f.seek(byte_range[0])
                if byte_range[1] is None:
                    bin = f.read()
                else:
                    bin = f.read(byte_range[1])

    ## file size
    fsize = len(bin)

    return(bin)
