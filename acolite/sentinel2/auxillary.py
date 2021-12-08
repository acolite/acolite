## def auxillary
## reads auxillary data included from S2 PB004 onward
## written by Quinten Vanhellemont, RBINS
## 2021-11-27
## modifications:

def auxillary(bundle, granule, sources = ['AUX_CAMSFO', 'AUX_ECMWFT'], key_name = 'cfVarName', reshape = False):
    import os, pygrib

    data = {}
    for source in sources:
        ## find aux grib files
        aux_file = '{}/GRANULE/{}/AUX_DATA/{}'.format(bundle, granule, source)
        if not os.path.exists(aux_file): continue

        ## open grib file
        with pygrib.open(aux_file) as gr:
            for ig, g in enumerate(gr):
                grb = gr.select()[ig]
                grb['stepRange'] = 's' # change stepRange to avoid errors - is "s" correct?

                ## read datasets for this parameter
                data[grb[key_name]] = {}
                for k in grb.keys():
                    try:
                        data[grb[key_name]][k] = grb[k]
                    except:
                        pass

                ## reshape longitudes and latitudes
                if reshape:
                    data[grb[key_name]]['longitudes'] = data[grb[key_name]]['longitudes'].reshape(int(grb['Ni']),int(grb['Nj']))
                    data[grb[key_name]]['latitudes'] = data[grb[key_name]]['latitudes'].reshape(int(grb['Ni']),int(grb['Nj']))
    return(data)
