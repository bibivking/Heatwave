#!/usr/bin/python

def leap_year(year):
   if (year % 4) == 0:
       return 366
   else:
       return 365
       
def get_reverse_colormap(var_name):

    '''
    To tell whether it needs a reversed colormap
    '''

    var_reverse_yes = [ "Rainf_f_inst","Rainf_tavg","Evap_tavg","ECanop_tavg","TVeg_tavg","ESoil_tavg",
                        "Qs_tavg","Qsb_tavg", "Snowf_tavg","GPP_tavg","Qle_tavg","SoilMoist_inst",
                        "FWsoil_tavg","SnowCover_inst","Qair_f_inst","Wind_f_inst","SWE_inst",
                        "SnowDepth_inst","SoilWet_inst", 
                        "rh2"]
    var_reverse_no  = [ "Qh_tavg","Qg_tavg","Swnet_tavg","Lwnet_tavg","SWdown_f_inst","LWdown_f_inst",
                        "VegT_tavg","AvgSurfT_tavg","Tair_f_inst","SoilTemp_inst","Albedo_inst",
                        "Psurf_f_inst",
                        "T2"]   

    if var_name in var_reverse_yes:
        return_value = True
    elif var_name in var_reverse_no:
        return_value = False
    else: 
        return_value = None

    return return_value