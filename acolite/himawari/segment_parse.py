## def segment read
## parses Himawari segment header and data
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-19
## modifications:

def segment_parse(segment_file, parse_header = True, parse_data = True, header = None):
    import struct
    import numpy as np
    import acolite as ac

    ## get total header length
    bin = ac.himawari.segment_read(segment_file, (70, 4))
    ret = ac.himawari.block_unpack(bin, [['total_header_length', 4, 'I']], strip = True)
    header_length = ret['total_header_length']

    ## get total data length
    bin = ac.himawari.segment_read(segment_file, (74, 4))
    ret = ac.himawari.block_unpack(bin, [['total_data_length', 4, 'I']], strip = True)
    data_length = ret['total_data_length']

    ## also parse header if parse_data is set
    if (header is not None): parse_header = False
    if (parse_data) & (header is None): parse_header = True

    if parse_header:
        ## read header bin data
        header_bin = ac.himawari.segment_read(segment_file, (0, header_length))
        ## parse header
        header = ac.himawari.header_parse(header_bin)

    if not parse_data: return(header)

    ## read data
    data_bin = ac.himawari.segment_read(segment_file, (header_length, data_length))

    ## get data dimensions from header
    dim = header[2]['number_of_columns'],header[2]['number_of_lines']
    data = struct.unpack('<{}H'.format(int(len(data_bin)/2)), data_bin)
    data = np.asarray(data).reshape(dim[1], dim[0])

    ## return data
    if not parse_header: return(data)
    return(header, data)
