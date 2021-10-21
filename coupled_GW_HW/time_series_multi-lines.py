#!/usr/bin/python

'''
Functions:
1. plot multi lines in one figure
'''

from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from matplotlib import cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_scale_offline
from common_utils import *

def plot_time_series( file_paths, var_name, date_s, date_e, loc_lat=None, loc_lon=None,
                      lat_name=None, lon_name=None, message=None,labal_names=None ):

    print("======== In plot_time_series =========")

    fig, ax = plt.subplots()
    # colors = cm.get_cmap('Paired')[:]
    colors = cm.Paired(np.arange(0,12))
    print(colors)
    # plot line 1
    
    delta = datetime(2009,1,28) - datetime(2000,1,1)

    for i, file_path in enumerate(file_paths):

        print(file_path)
        
        Time, Var = read_var(file_path, var_name, loc_lat, loc_lon, lat_name, lon_name)
        time, var = time_series_var(Time, Var, date_s, date_e)

        # t = []
        # for j in np.arange(len(time)):
        #     t.append(time[j].days+(time[j].seconds)/24./60./60. - delta.days)

        # ======== Local time 10am -4pm =======
        time_cood = []
        for j in np.arange(len(time)):
            time_cood.append( (time[j].seconds <= 6.*60.*60.) ) 

        days  = []
        value = []
        for j in np.arange(len(time)):
            if time_cood[j]:
                days.append(time[j].days)
                value.append(var[j])
        print(days)
        print(value)

        data          = pd.DataFrame(days, columns=['day'])
        data['value'] = value

        if var_name in ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Rainf_tavg","Qs_tavg","Qsb_tavg"]:
            data_group = data.groupby(by=["day"]).sum()
            scale = 1800.
            ax.plot(data_group.index-delta.days, data_group['value']*scale, c = colors[i], label=labal_names[i])
        elif var_name in ["AvgSurfT_tavg","VegT_tavg","Tair_f_inst"]:
            data_group = data.groupby(by=["day"]).mean()
            scale = -273.15
            ax.plot(data_group.index-delta.days, data_group['value']-scale, c = colors[i], label=labal_names[i])
        else:
            data_group = data.groupby(by=["day"]).mean()
            ax.plot(data_group.index-delta.days, data_group['value'], c = colors[i], label=labal_names[i])

        Time = None
        time = None
        Var  = None
        var  = None
        time_cood = None
        days = None
        value= None
        data = None
        data_group = None

    # ax.set_xlim([np.min(var1*scale,var2*scale), np.max(var1*scale,var2*scale)])
    # Time2, Var2 = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_name[1], lon_name[1])
    # time2, var2 = time_series_var(Time2,Var2,date_s,date_e)
    # var = np.zeros((2,len(var1)))
    # var[0,:] = var1
    # var[1,:] = var2
    # ax.plot(t1, var*scale, alpha=0.5)
    # ax.set_ylabel('mm')
    ax.set_title(var_name)
    # ax.set_xticks(x1[::1440])
    # ax.set_xticklabels(np.arange(date_s,date_e,1))
    ax.legend()

    fig.tight_layout()
    if message == None:
        message = var_name
    else:
        message = message + "_" + var_name
    if loc_lat != None:
        message = message + "_lat="+str(loc_lat) + "_lon="+str(loc_lon)

    plt.savefig('./plots/19Oct/time_series_multi-lines_'+message+'.png',dpi=300)

if __name__ == "__main__":

    # #######################
    #        Variables      #
    # #######################

    var_names  = [  "Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg",
                    "Qh_tavg","Qle_tavg","Qg_tavg","Qs_tavg","Qsb_tavg",
                    "FWsoil_tavg","AvgSurfT_tavg","VegT_tavg","Tair_f_inst",
                    "Rainf_tavg", "Qair_f_inst" ]

    # #######################
    #     path setting      #
    # #######################

    hw_name   = "hw2009_15Oct"
            #  [  "hw2009_15Oct", "hw2011_15Oct", "hw2013_15Oct", "hw2014_15Oct","hw2017_15Oct", "hw2019_15Oct" ] 
    path      = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/" + hw_name
    rst_dates = ["20090113",
                 "20090115",
                 "20090117",
                 "20090119",
                 "20090121" ]
    periods   = ["20090122-20090213"]

    file_paths = [  path + "/fd_rst_" + rst_dates[0] + "/LIS_output/LIS.CABLE." + periods[0] + ".nc",
                    path + "/gw_rst_" + rst_dates[0] + "/LIS_output/LIS.CABLE." + periods[0] + ".nc",
                    path + "/fd_rst_" + rst_dates[1] + "/LIS_output/LIS.CABLE." + periods[0] + ".nc",
                    path + "/gw_rst_" + rst_dates[1] + "/LIS_output/LIS.CABLE." + periods[0] + ".nc",
                    path + "/fd_rst_" + rst_dates[2] + "/LIS_output/LIS.CABLE." + periods[0] + ".nc",
                    path + "/gw_rst_" + rst_dates[2] + "/LIS_output/LIS.CABLE." + periods[0] + ".nc",
                    path + "/fd_rst_" + rst_dates[3] + "/LIS_output/LIS.CABLE." + periods[0] + ".nc",
                    path + "/gw_rst_" + rst_dates[3] + "/LIS_output/LIS.CABLE." + periods[0] + ".nc",
                    path + "/fd_rst_" + rst_dates[4] + "/LIS_output/LIS.CABLE." + periods[0] + ".nc",
                    path + "/gw_rst_" + rst_dates[4] + "/LIS_output/LIS.CABLE." + periods[0] + ".nc" ] 

 
    wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/fd_feb2009/WRF_output/wrfout_d01_2009-02-01_00:00:00"

    ####################################
    #         plot_time_series         #
    ####################################

    date_s       = datetime(2009,1,22,0,0,0)
    date_e       = datetime(2009,2,13,23,59,59)
    loc_lat      = [-35.25,-34.75]
    loc_lon      = [149.75,150.25]
    lat_name     = "lat"
    lon_name     = "lon"
    labal_names  = [ "fd-14d","gw-14d","fd-12d","gw-12d","fd-10d",
                     "gw-10d","fd-8d","gw-8d","fd-6d","gw-6d" ]
    for var_name in var_names:
        message   = "hw2009_15Oct"
        plot_time_series(file_paths, var_name, date_s, date_e, loc_lat=loc_lat, loc_lon=loc_lon,
                         lat_name=lat_name, lon_name=lon_name, message=message,labal_names=labal_names)

