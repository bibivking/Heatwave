#!/usr/bin/python

def get_var_scale(var_name):

    '''
    Convert units
    '''

    var_s2d        = ["Rainf_f_inst","Rainf_tavg","Evap_tavg","ECanop_tavg","TVeg_tavg","ESoil_tavg","Qs_tavg","Qsb_tavg",
                      "Snowf_tavg"]
    var_umol_s2g_d = ["GPP_tavg"]
    var_wm2        = ["Qle_tavg","Qh_tavg","Qg_tavg","Swnet_tavg","Lwnet_tavg","SWdown_f_inst","LWdown_f_inst"]
    var_degc       = ["VegT_tavg","AvgSurfT_tavg","Tair_f_inst","SoilTemp_inst"]
    var_percent    = ["Albedo_inst","FWsoil_tavg","SnowCover_inst","Qair_f_inst","SoilMoist_inst"]
    var_ms         = ["Wind_f_inst"]
    var_hPa        = ["Psurf_f_inst"]
    var_mm         = ["SWE_inst"]
    var_m          = ["SnowDepth_inst","SoilWet_inst"]

    s2d            = 3600*24. # s-1 to d-1
    mmd2wm2        = 28.94    # 1 mm/day = 28.94 W/m2
    umol_s2g_d     = 0.000001*12*s2d # umol s-1 to g d-1
    Pa2hPa         = 0.01     # 1 Pa = 0.01 hPa
    degK2C         = -273.15  # 1 K = -273.15 C

    if var_name in var_s2d:
        scale = s2d*mmd2wm2
        units = "W/m2" #"mm d-1"
    elif var_name in var_umol_s2g_d:
        scale = umol_s2g_d
        units = "gC m-2 d-1"
    elif var_name in var_wm2:
        scale = 1.
        units = "W/m2"
    elif var_name in var_degc:
        scale = degK2C
        units = "degC"
    elif var_name in var_percent:
        scale = 1.
        units = "-"
    elif var_name in var_ms:
        scale = 1.
        units = "m/s"
    elif var_name in var_hPa:
        scale = Pa2hPa
        units = "hPa"
    elif var_name in var_mm:
        scale = 1.
        units = "mm"
    elif var_name in var_m:
        scale = 1.
        units = "m"
    else:
        scale = 1.
        units = None
    return (scale, units)

def get_var_range_diff(var_name):

    '''
    Convert units
    '''

    var_s2d        = ["Rainf_f_inst","Rainf_tavg","Evap_tavg","ECanop_tavg","TVeg_tavg","ESoil_tavg","Qs_tavg","Qsb_tavg",
                      "Snowf_tavg"]
    var_umol_s2g_d = ["GPP_tavg"]
    var_wm2        = ["Qle_tavg","Qh_tavg","Qg_tavg","Swnet_tavg","Lwnet_tavg","SWdown_f_inst","LWdown_f_inst"]
    var_degc       = ["VegT_tavg","AvgSurfT_tavg","Tair_f_inst","SoilTemp_inst"]
    var_percent    = ["Albedo_inst","FWsoil_tavg","SnowCover_inst","Qair_f_inst","SoilMoist_inst"]
    var_ms         = ["Wind_f_inst"]
    var_hPa        = ["Psurf_f_inst"]
    var_mm         = ["SWE_inst"]
    var_m          = ["SnowDepth_inst","SoilWet_inst"]

    ranges         = [0.0,0.0]

    if var_name in var_s2d:
        ranges[0] = -100.
        ranges[1] = 100.
    elif var_name in var_umol_s2g_d:
        ranges[0] = -20.
        ranges[1] = 20.
    elif var_name in var_wm2:
        ranges[0] = -100.
        ranges[1] = 100.
    elif var_name in var_degc:
        ranges[0] = -5.
        ranges[1] = 5.
    elif var_name in var_percent:
        ranges[0] = -0.5
        ranges[1] = 0.5
    elif var_name in var_ms:
        ranges[0] = -10.
        ranges[1] = 10.
    elif var_name in var_hPa:
        ranges[0] = -50.
        ranges[1] = 50.
    elif var_name in var_mm:
        ranges[0] = -100.
        ranges[1] = 100.
    elif var_name in var_m:
        ranges[0] = -1.
        ranges[1] = 1.
    else:
        ranges = None
    return ranges