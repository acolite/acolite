## ancillary_list
## lists ancillary data from a given date from the ocean data server
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-10-17
## modifications: 2021-01-31 (QV) added 2020-05-28 limit for TOAST
##                2021-03-01 (QV) simplified for acg renamed from ancillary_list

def list_files(date):
    import os
    from datetime import datetime, timedelta

    year, month, day = [int(i) for i in date.split('-')]
    dtime = datetime(year, month, day)

    isodate = dtime.strftime("%Y-%m-%d")
    year = dtime.strftime("%Y")
    jday = dtime.strftime("%j").zfill(3)
    yjd = "{}{}".format(dtime.strftime("%Y"),jday)

    dtime_next = dtime+timedelta(days=1)
    year_next = dtime_next.strftime("%Y")
    jday_next = dtime_next.strftime("%j").zfill(3)
    yjd_next = "{}{}".format(year_next,jday_next)

    if isodate > '2020-05-27':
        file_types = ['MET_NCEPR2_6h','MET_NCEPR2_6h_NEXT']
    else:
        file_types = ['TOAST','MET_NCEPR2_6h','MET_NCEPR2_6h_NEXT']

    if ((int(year) == 2004) & (int(jday) > 336)) | (int(year) > 2004):
        file_types+=['O3_AURAOMI_24h']
    else:
        file_types+=['O3_EPTOMS_24h', 'O3_N7TOMS_24h']

    basefiles = []
    for file_type in file_types:
        if file_type == "TOAST":
            cfile = "S{}00{}23_TOAST.OZONE".format(yjd,jday)
            basefiles.append(cfile)
        elif file_type == "O3_AURAOMI_24h":
            cfile = "N{}00_O3_AURAOMI_24h.hdf".format(yjd)
            basefiles.append(cfile)
        elif file_type == "O3_EPTOMS_24h":
            cfile = "N{}00_O3_EPTOMS_24h.hdf".format(yjd)
            basefiles.append(cfile)
        elif file_type == "O3_N7TOMS_24h":
            cfile = "N{}00_O3_N7TOMS_24h.hdf".format(yjd)
            basefiles.append(cfile)
        elif file_type == "MET_NCEP":
            cfile = ["S{}{}_NCEP.MET".format(yjd,h) for h in ['00','06','12','18']]
            for f in cfile: basefiles.append(f)
        elif file_type == "MET_NCEP_NEXT":
            cfile = "S{}00_NCEP.MET".format(yjd_next)
            basefiles.append(cfile)
        elif file_type == "MET_NCEP_6h":
            cfile = ["N{}{}_MET_NCEP_6h.hdf".format(yjd,h) for h in ['00','06','12','18']]
            for f in cfile: basefiles.append(f)
        elif file_type == "MET_NCEP_6h_NEXT":
            cfile = "N{}00_MET_NCEP_6h.hdf".format(yjd_next)
            basefiles.append(cfile)
        elif file_type == "MET_NCEPR2_6h":
            cfile = ["N{}{}_MET_NCEPR2_6h.hdf.bz2".format(yjd,h) for h in ['00','06','12','18']]
            for f in cfile: basefiles.append(f)
        elif file_type == "MET_NCEPR2_6h_NEXT":
            cfile = "N{}00_MET_NCEPR2_6h.hdf.bz2".format(yjd_next)
            basefiles.append(cfile)
        else:
            continue
    return(basefiles)
