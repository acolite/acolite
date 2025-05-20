## def block_unpack
## unpacks Himawari HSD data block
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-19
## modifications:

def block_unpack(bin, block, strip = True, encoding = 'ascii'):
    import struct

    ret = {}

    ## unpack binary data
    pos = 0
    for v in block:
        d = bin[pos:pos+v[1]]
        pos += v[1]
        if v[2] == '%ds':
            res = struct.unpack(v[2] % v[1], d)[0].decode(encoding)
            if strip: res = res.strip('\x00')
        else:
            res = (struct.unpack("<{}".format(v[2]), d)[0])
        ret[v[0]] = res

    return(ret)
