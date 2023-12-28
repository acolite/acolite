## ancillary_list
## lists ancillary data from a given date from the ocean data server
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-10-17
## modifications: 2021-01-31 (QV) added 2020-05-28 limit for TOAST
##                2021-03-01 (QV) simplified for acg renamed from ancillary_list
##                2023-10-16 (QV) changed date parsing, added GMAO FP and MERRA2 files
##                2023-12-28 (QV) added GMAO_IT_MET

def list_files(date, file_types = None):
    import os
    import acolite as ac
    import datetime, dateutil.parser
    dtime = dateutil.parser.parse(date)

    isodate = dtime.strftime("%Y-%m-%d")
    year = dtime.strftime("%Y")
    month = dtime.strftime("%m")
    day = dtime.strftime("%d")

    jday = dtime.strftime("%j").zfill(3)
    yjd = "{}{}".format(dtime.strftime("%Y"),jday)

    dtime_next = dtime+datetime.timedelta(days=1)
    year_next = dtime_next.strftime("%Y")
    jday_next = dtime_next.strftime("%j").zfill(3)
    yjd_next = "{}{}".format(year_next,jday_next)

    ## get ancillary file types
    if file_types is None:
        file_types = ac.settings['run']['ancillary_type']
    if type(file_types) is not list: file_types = [file_types]

    # if isodate > '2020-05-27':
    #     file_types = ['MET_NCEPR2_6h','MET_NCEPR2_6h_NEXT']
    # else:
    #     file_types = ['TOAST','MET_NCEPR2_6h','MET_NCEPR2_6h_NEXT']
    #
    # file_types += ['MET_NCEP_6h','MET_NCEP_6h_NEXT']
    # file_types += ['GMAO_MERRA2_MET']
    #
    # if ((int(year) == 2004) & (int(jday) > 336)) | (int(year) > 2004):
    #     file_types+=['O3_AURAOMI_24h']
    # else:
    #     file_types+=['O3_EPTOMS_24h', 'O3_N7TOMS_24h']

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
        elif file_type == "GMAO_MERRA2_MET":
            ## hourly file before
            cfile = "GMAO_MERRA2.{}T{}0000.MET.nc".format(dtime.strftime("%Y%m%d"),str(dtime.hour).zfill(2))
            basefiles.append(cfile)
            ## hourly file after
            if dtime.hour < 23:
                cfile = "GMAO_MERRA2.{}T{}0000.MET.nc".format(dtime.strftime("%Y%m%d"),str(dtime.hour+1).zfill(2))
            else:
                cfile = "GMAO_MERRA2.{}T{}0000.MET.nc".format((dtime+datetime.timedelta(days=1)).strftime("%Y%m%d"),'00')
            basefiles.append(cfile)
        elif file_type == "GMAO_FP_MET":
            ## 3 hourly file before
            h0 = dtime.hour - (dtime.hour % 3)
            cfile = "GMAO_FP.{}T{}0000.MET.NRT.nc".format(dtime.strftime("%Y%m%d"),str(h0).zfill(2))
            basefiles.append(cfile)
            ## 3 hourly file after
            if h0 != 21:
                cfile = "GMAO_FP.{}T{}0000.MET.NRT.nc".format(dtime.strftime("%Y%m%d"),str(h0+3).zfill(2))
            else:
                cfile = "GMAO_FP.{}T{}0000.MET.NRT.nc".format((dtime+datetime.timedelta(days=1)).strftime("%Y%m%d"),'00')
            basefiles.append(cfile)
        elif file_type == "GMAO_IT_MET":
            ## hourly file before
            cfile = "GMAO_IT.{}T{}0000.MET.NRT.nc".format(dtime.strftime("%Y%m%d"),str(dtime.hour).zfill(2))
            basefiles.append(cfile)
            ## hourly file after
            if dtime.hour < 23:
                cfile = "GMAO_IT.{}T{}0000.MET.NRT.nc".format(dtime.strftime("%Y%m%d"),str(dtime.hour+1).zfill(2))
            else:
                cfile = "GMAO_IT.{}T{}0000.MET.NRT.nc".format((dtime+datetime.timedelta(days=1)).strftime("%Y%m%d"),'00')
            basefiles.append(cfile)
        else:
            continue
    return(basefiles)
