## def gpt_geometry
## computes per pixel geometry using SNAP gpt
##
## written by Quinten Vanhellemont, RBINS
## 2021-02-17
## modifications: 2021-02-23 (QV) added quotation marks to gpt command to support spaces in path names

def gpt_geometry(bundle, output=None, target_res=60, override=True, verbosity=0, format='GeoTIFF'):
    import os
    import acolite as ac
    import subprocess
    files = []
    res = '{:.0f}'.format(target_res)

    gpt = '{}/bin/gpt'.format(ac.config['snap_directory'])
    if os.path.exists(gpt) is False:
        if verbosity>0: print('gpt not found at {}'.format(gpt))
        return(files)

    if output is None: output = os.path.dirname(bundle)
    if not os.path.exists(output): os.makedirs(output)

    parameters = ['view_zenith_mean','view_azimuth_mean','sun_zenith','sun_azimuth']
    if format == 'GeoTIFF':
        ext = 'tif'
    elif format == 'NetCDF4-BEAM':
        ext = 'nc'
        parameters=[','.join(parameters)]
    elif format == 'NetCDF4-CF':
        ext = 'nc'
        parameters=[','.join(parameters)]
    else:
        if verbosity>0: print('format {} not configures'.format(format))
        return(files)

    for parameter in parameters:
        if format == 'GeoTIFF':
            geometry_file = '{}/{}'.format(output, os.path.basename(bundle).replace('.SAFE', '_geometry_{}m_{}.{}'.format(res, parameter, ext)))
            gptfile = '{}/gpt_geometry_graph_{}m_{}.xml'.format(output,res,parameter)
        elif format == 'NetCDF4-BEAM':
            geometry_file = '{}/{}'.format(output, os.path.basename(bundle).replace('.SAFE', '_geometry_{}m.{}'.format(res, ext)))
            gptfile = '{}/gpt_geometry_graph_{}m.xml'.format(output,res)
        elif format == 'NetCDF4-CF':
            geometry_file = '{}/{}'.format(output, os.path.basename(bundle).replace('.SAFE', '_geometry_{}m.{}'.format(res, ext)))
            gptfile = '{}/gpt_geometry_graph_{}m.xml'.format(output,res)

        if not os.path.exists(geometry_file) or override:
            ## create gpt graph
            ifile = '{}/S2/gpt_geometry_graph.xml'.format(ac.config['data_dir'])
            if verbosity > 0: print('Writing gpt graph to {}'.format(gptfile, res))
            with open(ifile, 'r') as fi, open(gptfile, 'w') as fo:
                for line in fi.readlines():
                    if '$bundle' in line:
                        line = line.replace('$bundle', bundle)
                    if '$target_res' in line:
                        line = line.replace('$target_res', res)
                    if '$geometry_file' in line:
                        line = line.replace('$geometry_file', geometry_file)
                    if '$parameter' in line:
                        line = line.replace('$parameter', parameter)
                    if '$format' in line:
                        line = line.replace('$format', format)
                    fo.write(line)
            ## run processing
            if verbosity > 0: print('Running gpt resampling to {}m'.format(res))
            sp = subprocess.run("'{}'".format(gpt) + " " + "'{}'".format(gptfile), shell=True, check=True, stdout=subprocess.PIPE)

        if os.path.exists(geometry_file): files.append(geometry_file)
    return(files)
