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
from matplotlib.cm import get_cmap
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_scale_offline
from common_utils import *

def plot_time_series( file_paths, var_name, date_s, date_e, hw_start, loc_lat=None, loc_lon=None,
                      lat_name=None, lon_name=None, message=None,labal_names=None ):

    print("======== In plot_time_series =========")

    fig, ax = plt.subplots()
    # colors = cm.get_cmap('Paired')[:]
    colors = cm.Paired(np.arange(0,12))
    print(colors)
    # plot line 1

    delta = hw_start - datetime(2000,1,1)

    for i, file_path in enumerate(file_paths):

        print(file_path)

        Time, Var = read_var(file_path, var_name, loc_lat, loc_lon, lat_name, lon_name)
        time, var = time_series_var(Time, Var, date_s, date_e)
        # print(time[:2])
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
            scale = 273.15
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

    plt.savefig('./plots/5Nov/hw_evolution/3Nov/time_series_multi-lines_'+message+'.png',dpi=300)

def plot_time_series_errorbar( file_paths, var_name, date_s, date_e, hw_start, loc_lat=None, loc_lon=None,
                      lat_name=None, lon_name=None, message=None,labal_names=None ):

    print("======== In plot_time_series_errorbar =========")

    fig, ax = plt.subplots()
    colors = cm.Paired(np.arange(0,12))
    # print(colors)

    delta = hw_start - datetime(2000,1,1)

    total_days = (date_e-date_s).days + 1
    # print(total_days)

    values = np.zeros((len(file_paths),total_days))

    for i, file_path in enumerate(file_paths):

        # print(file_path)

        Time, Var = read_var(file_path, var_name, loc_lat, loc_lon, lat_name, lon_name)
        time, var = time_series_var(Time, Var, date_s, date_e)

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
        # print(days)
        # print(value)

        data          = pd.DataFrame(days, columns=['day'])
        data['value'] = value

        if var_name in ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Rainf_tavg","Qs_tavg","Qsb_tavg"]:
            data_group = data.groupby(by=["day"]).sum()
            scale = 1800.
            values[i,:] = data_group['value']*scale
        elif var_name in ["AvgSurfT_tavg","VegT_tavg","Tair_f_inst"]:
            data_group = data.groupby(by=["day"]).mean()
            scale = 273.15
            values[i,:] = data_group['value'] - scale
        else:
            data_group = data.groupby(by=["day"]).mean()
            values[i,:] = data_group['value']
        if i == 0:
            times = data_group.index-delta.days

        Time = None
        time = None
        Var  = None
        var  = None
        time_cood = None
        days = None
        value= None
        data = None
        data_group = None
        rain_mask_multi = None
        wtd_mask_multi  = None
        pft_mask_multi  = None

    gw_mean = np.mean(values[0:5,:], axis =0)
    gw_min  = np.min(values[0:5,:], axis =0)
    gw_max  = np.max(values[0:5,:], axis =0)
    fd_mean = np.mean(values[5:10,:], axis =0)
    fd_min  = np.min(values[5:10,:], axis =0)
    fd_max  = np.max(values[5:10,:], axis =0)

    ax.plot(times, gw_mean, color="red", ls='-', lw=2.0, zorder=10, label="FD")
    ax.fill_between(times, gw_min, gw_max, alpha=0.5, facecolor='red', zorder=10)

    ax.plot(times, fd_mean, color="blue", ls='-', lw=2.0, zorder=10, label="GW")
    ax.fill_between(times, fd_min, fd_max, alpha=0.5, facecolor='blue', zorder=10)

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

    plt.savefig('./plots/19Oct/hw_evolution/time_series_multi-lines_errorbar_'+message+'.png',dpi=300)

def plot_time_series_errorbar_select_regions( file_paths, var_name, date_s, date_e, seconds, hw_start, loc_lat=None, loc_lon=None,
                      lat_name=None, lon_name=None, message=None,labal_names=None, rain_val=None,wtd_val=None,
                      pft_val=None, ensemble_path=None, check_pixel=False ):

    print("======== In plot_time_series_errorbar =========")

    fig, ax = plt.subplots()
    colors = cm.Paired(np.arange(0,12))

    delta = hw_start - datetime(2000,1,1)

    total_days = (date_e-date_s).days + 1

    values_mean = np.zeros((len(file_paths),total_days))
    values_max  = np.zeros((len(file_paths),total_days))
    values_min  = np.zeros((len(file_paths),total_days))

    if rain_val is not None:
        Time, rain= read_var(ensemble_path, "Rainf_tavg", loc_lat, loc_lon, lat_name, lon_name)
        rain_mean = spital_var(UTC_to_AEST(Time), rain, date_s, date_e)
        rain_mask = np.ones(np.shape(rain_mean), dtype=bool)
        rain_mask = np.where( rain_mean*24.*3600. < rain_val, True, False)
        print("(rain_mask == True).sum()")
        print((rain_mask == True).sum())

    if wtd_val is not None:
        Time, wtd = read_var(ensemble_path, "WaterTableD_tavg", loc_lat, loc_lon, lat_name, lon_name)
        wtd_mean = spital_var(UTC_to_AEST(Time), wtd, date_s, date_e)
        wtd_mask = np.ones(np.shape(wtd_mean), dtype=bool)
        wtd_mask = np.where(np.all([(wtd_mean/1000.) >=wtd_val[0], (wtd_mean/1000.) <wtd_val[1]],axis=0), True, False)
        print("(wtd_mask == True).sum()")
        print((wtd_mask == True).sum())

    if pft_val is not None:
        pft      = tree_mask(ensemble_path,"Landcover_inst")
        pft_mask = np.ones(np.shape(pft), dtype=bool)
        if len(pft_val) == 1:
            pft_mask = np.where( pft == pft_val[0], True, False)
        elif (len(pft_val) == 2) :
            pft_mask = np.where( np.all([pft >=pft_val[0], pft<=pft_val[1]],axis=0), True, False)

    var = []
    for i, file_path in enumerate(file_paths):

        Time, Var = read_var(file_path, var_name, loc_lat, loc_lon, lat_name, lon_name)
        ntime = len(Time)

        if rain_mask is not None:
            rain_mask_multi = [rain_mask]*ntime
            # Var = np.where(rain_mask_multi == True, Var, np.nan)
        if wtd_mask is not None:
            wtd_mask_multi  = [wtd_mask]*ntime
            # Var = np.where(wtd_mask_multi == True, Var, np.nan)
        if pft_mask is not None:
            pft_mask_multi  = [pft_mask]*ntime
            # Var = np.where(pft_mask_multi == True, Var, np.nan)

        Var_tmp = np.where( np.all([rain_mask_multi,wtd_mask_multi,pft_mask_multi],axis=0), Var, np.nan)
        # print(np.any(np.isnan(Var_tmp[0,:,:])))


        # check pixel cover
        # if check_pixel:
        #
        #     # check mask
        #     file_tmp = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg/wrfout_20090122-20090213_gw"
        #     file = nc.Dataset(file_tmp, mode='r')
        #     p1   = getvar(file, "pressure")
        #     lon  = file.variables['XLONG'][0,:,:]
        #     lat  = file.variables['XLAT'][0,:,:]
        #
        #     fig1 = plt.figure(figsize=(12,9))
        #     # Get the lat/lon coordinates
        #
        #     ax1 = plt.axes(projection=ccrs.PlateCarree())
        #     ax1.set_extent([110,155,-45,-10])
        #     ax1.coastlines(resolution="50m",linewidth=1)
        #
        #     # Add gridlines
        #     gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        #     gl.xlabels_top  = False
        #     gl.ylabels_right= False
        #     gl.xlines       = True
        #     gl.xlocator     = mticker.FixedLocator(np.arange(110,155,1))
        #     gl.ylocator     = mticker.FixedLocator(np.arange(-45,-10,1))
        #     gl.xformatter   = LONGITUDE_FORMATTER
        #     gl.yformatter   = LATITUDE_FORMATTER
        #     gl.xlabel_style = {'size':10, 'color':'black','rotation': 90}
        #     gl.ylabel_style = {'size':10, 'color':'black'}
        #
        #     cmap = plt.cm.jet
        #
        #     mask = np.where(np.isnan(np.nanmean(Var_tmp,axis=0)), 0, 1)
        #     plt.contourf(lon, lat, mask, levels=[-0.1,0.1,1.1], transform=ccrs.PlateCarree(), cmap=cmap) #,"jet" #“rainbow”#"coolwarm"
        #     # plt.pcolormesh( mask, transform=ccrs.PlateCarree(),cmap=cmap)       # lon, lat,levels=[-0.1,0.1,1.1],
        #     plt.savefig('./plots/5Nov/pixel_selection/check_pixels_'+message+'.png',dpi=300)
        #     # np.savetxt(message,mask)

        var_mean, var_max, var_min = time_series_statistic(UTC_to_AEST(Time), Var_tmp, date_s, date_e) # shape = (selected_pixels_num, time_series)

        # ======== Local time 10am -4pm =======
        time = time_series_time(UTC_to_AEST(Time), date_s, date_e)

        if seconds == None:
            time_cood = time
        else:
            time_cood = []
            for j in np.arange(len(time)):
                if seconds[0] > seconds[1]:
                    time_cood.append( (time[j].seconds >= seconds[0]) | (time[j].seconds < seconds[1]) )
                elif seconds[0] < seconds[1]:
                    time_cood.append( (time[j].seconds >= seconds[0]) & (time[j].seconds < seconds[1]) )
        # print(time_cood)

        days       = []
        value_mean = []
        value_max  = []
        value_min  = []
        for j in np.arange(len(time)):
            if time_cood[j]:
                days.append(time[j].days)
                value_mean.append(var_mean[j])
                value_max.append(var_max[j])
                value_min.append(var_min[j])

        data          = pd.DataFrame(days, columns=['day'])
        data['mean']  = value_mean
        data['max']   = value_max
        data['min']   = value_min

        if var_name in ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Rainf_tavg","Qs_tavg","Qsb_tavg"]:
            data_group = data.groupby(by=["day"]).sum()
            scale = 1800.
            values_mean[i,:] = data_group['mean']*scale
            values_max[i,:]  = data_group['max']*scale
            values_min[i,:]  = data_group['min']*scale
        elif var_name in ["AvgSurfT_tavg","VegT_tavg","Tair_f_inst"]:
            data_group = data.groupby(by=["day"]).mean()
            scale = 273.15
            values_mean[i,:] = data_group['mean'] - scale
            values_max[i,:]  = data_group['max'] - scale
            values_min[i,:]  = data_group['min'] - scale
        else:
            data_group = data.groupby(by=["day"]).mean()
            values_mean[i,:] = data_group['mean']
            values_max[i,:]  = data_group['max']
            values_min[i,:]  = data_group['min']

        if i == 0:
            times = data_group.index-delta.days

        Time = None
        time = None
        Var  = None
        var  = None
        time_cood = None
        days = None
        value= None
        data = None
        data_group = None
        rain_mask_multi = None
        wtd_mask_multi  = None
        pft_mask_multi  = None

    gw_mean = np.mean(values_mean[0:5,:], axis =0)
    gw_max  = np.max(values_max[0:5,:], axis =0)
    gw_min  = np.min(values_min[0:5,:], axis =0)
    fd_mean = np.mean(values_mean[5:10,:], axis =0)
    fd_max  = np.max(values_max[5:10,:], axis =0)
    fd_min  = np.min(values_min[5:10,:], axis =0)

    ax.plot(times, gw_mean, color="red", ls='-', lw=2.0, zorder=10, label="FD")
    ax.fill_between(times, gw_min, gw_max, alpha=0.5, facecolor='red', zorder=10)

    ax.plot(times, fd_mean, color="blue", ls='-', lw=2.0, zorder=10, label="GW")
    ax.fill_between(times, fd_min, fd_max, alpha=0.5, facecolor='blue', zorder=10)

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

    plt.savefig('./plots/5Nov/hw_evolution/3Nov/time_series_multi-lines_'+message+'.png',dpi=300)

if __name__ == "__main__":

    # #######################
    #        Variables      #
    # #######################
    #
    # var_names  = [  "Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg",
    #                 "Qh_tavg","Qle_tavg","Qg_tavg","Qs_tavg","Qsb_tavg",
    #                 "FWsoil_tavg","AvgSurfT_tavg","VegT_tavg","Tair_f_inst",
    #                 "Rainf_tavg", "Qair_f_inst" ]

    var_names  = [ "Evap_tavg", "TVeg_tavg","Qh_tavg","Qle_tavg",
                   "FWsoil_tavg","AvgSurfT_tavg","VegT_tavg","Tair_f_inst",
                   "Rainf_tavg", "Qair_f_inst" ]
    # #
    # #######################
    #     path setting      #
    # #######################

    hw_name = "hw2013_3Nov" #"hw2009_3Nov"
    check_pixel = True

    # 2009
    if hw_name == "hw2009_3Nov":
        start_date= "20090122"
        rst_dates = ["20090117","20090118","20090119","20090120","20090121" ]
        hw_start  = datetime(2009,1,28,0,0,0)
        if check_pixel:
            end_date  = "20090213"
            date_s    = datetime(2009,1,22,0,0,0) #datetime(2009,1,17,0,0,0)
            date_e    = datetime(2009,2,13,23,59,0)
        else:
            end_date  = "20090214"
            date_s    = datetime(2009,1,17,0,0,0)
            date_e    = datetime(2009,2,13,23,59,0)

    elif hw_name == "hw2013_3Nov":
        start_date= "20121229"
        rst_dates = ["20121224","20121225","20121226","20121227","20121228" ]
        hw_start  = datetime(2013,1,4,0,0,0)
        if check_pixel:
            end_date  = "20130122"
            date_s    = datetime(2012,12,29,0,0,0) #datetime(2012,12,24,0,0,0)
            date_e    = datetime(2013,1,22,23,59,0)
        else:
            end_date  = "20130123"
            date_s    = datetime(2012,12,24,0,0,0)
            date_e    = datetime(2013,1,22,23,59,0)

    elif hw_name == "hw2019_3Nov":
        start_date= "20190108"
        rst_dates = ["20190103","20190104","20190105","20190106","20190107" ]
        hw_start  = datetime(2019,1,14,0,0,0)
        if check_pixel:
            end_date  = "20190130"
            date_s    = datetime(2019,1,8,0,0,0) #datetime(2019,1,3,0,0,0)
            date_e    = datetime(2019,1,30,23,59,0)
        else:
            end_date  = "20190131"
            date_s    = datetime(2019,1,3,0,0,0)
            date_e    = datetime(2019,1,30,23,59,0)

    # #######################
    #     path setting      #
    # #######################

    path          = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/" + hw_name
    ensemble_path = path + "/ensemble_avg/LIS.CABLE." + start_date + "-" + end_date+ "_gw.nc"
    if check_pixel:
        file_paths = [  path + "/fd_rst_" + rst_dates[0] + "/LIS_output/LIS.CABLE." + start_date + "-" + end_date+ ".nc",
                        path + "/fd_rst_" + rst_dates[1] + "/LIS_output/LIS.CABLE." + start_date + "-" + end_date+ ".nc",
                        path + "/fd_rst_" + rst_dates[2] + "/LIS_output/LIS.CABLE." + start_date + "-" + end_date+ ".nc",
                        path + "/fd_rst_" + rst_dates[3] + "/LIS_output/LIS.CABLE." + start_date + "-" + end_date+ ".nc",
                        path + "/fd_rst_" + rst_dates[4] + "/LIS_output/LIS.CABLE." + start_date + "-" + end_date+ ".nc",
                        path + "/gw_rst_" + rst_dates[0] + "/LIS_output/LIS.CABLE." + start_date + "-" + end_date+ ".nc",
                        path + "/gw_rst_" + rst_dates[1] + "/LIS_output/LIS.CABLE." + start_date + "-" + end_date+ ".nc",
                        path + "/gw_rst_" + rst_dates[2] + "/LIS_output/LIS.CABLE." + start_date + "-" + end_date+ ".nc",
                        path + "/gw_rst_" + rst_dates[3] + "/LIS_output/LIS.CABLE." + start_date + "-" + end_date+ ".nc",
                        path + "/gw_rst_" + rst_dates[4] + "/LIS_output/LIS.CABLE." + start_date + "-" + end_date+ ".nc" ]
    else:
        file_paths = [  path + "/fd_rst_" + rst_dates[0] + "/LIS_output/LIS.CABLE." + rst_dates[0] + "-" + end_date+ ".nc",
                        path + "/fd_rst_" + rst_dates[1] + "/LIS_output/LIS.CABLE." + rst_dates[1] + "-" + end_date+ ".nc",
                        path + "/fd_rst_" + rst_dates[2] + "/LIS_output/LIS.CABLE." + rst_dates[2] + "-" + end_date+ ".nc",
                        path + "/fd_rst_" + rst_dates[3] + "/LIS_output/LIS.CABLE." + rst_dates[3] + "-" + end_date+ ".nc",
                        path + "/fd_rst_" + rst_dates[4] + "/LIS_output/LIS.CABLE." + rst_dates[4] + "-" + end_date+ ".nc",
                        path + "/gw_rst_" + rst_dates[0] + "/LIS_output/LIS.CABLE." + rst_dates[0] + "-" + end_date+ ".nc",
                        path + "/gw_rst_" + rst_dates[1] + "/LIS_output/LIS.CABLE." + rst_dates[1] + "-" + end_date+ ".nc",
                        path + "/gw_rst_" + rst_dates[2] + "/LIS_output/LIS.CABLE." + rst_dates[2] + "-" + end_date+ ".nc",
                        path + "/gw_rst_" + rst_dates[3] + "/LIS_output/LIS.CABLE." + rst_dates[3] + "-" + end_date+ ".nc",
                        path + "/gw_rst_" + rst_dates[4] + "/LIS_output/LIS.CABLE." + rst_dates[4] + "-" + end_date+ ".nc" ]

    wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+hw_name+"/ensemble_avg/wrfout_20090122-20090213_gw"


    ####################################
    #         plot_time_series         #
    ####################################
    lat_name     = "lat"
    lon_name     = "lon"
    labal_names  = [ "fd-10d","gw-10d","fd-9d","gw-9d","fd-8d",
                     "gw-8d","fd-7d","gw-7d","fd-6d","gw-6d" ]
    # labal_names  = [ "fd-14d","gw-14d","fd-12d","gw-12d","fd-10d",
    #                  "gw-10d","fd-8d","gw-8d","fd-6d","gw-6d" ]

    seconds       = [6.*60.*60.,18.*60.*60.]
    if seconds[0] < seconds[1]:
        time_spell = "daytime"
    for var_name in var_names:

        loc_lat   = [-27,-26]
        loc_lon   = [144.,145.]

        message   = hw_name +"_grass_wtd=0-5_"+time_spell

        plot_time_series_errorbar_select_regions( file_paths, var_name, date_s, date_e, seconds, hw_start,
                    loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name, message=message,
                    labal_names=labal_names, rain_val=10., wtd_val=[0.,5.], pft_val=[6],
                    ensemble_path=ensemble_path,check_pixel=check_pixel )

        loc_lat   = [-33,-32]
        loc_lon   = [143.,144.]

        message   = hw_name +"_grass_wtd=5-10_"+time_spell

        plot_time_series_errorbar_select_regions( file_paths, var_name, date_s, date_e, seconds, hw_start,
                    loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name, message=message,
                    labal_names=labal_names, rain_val=10., wtd_val=[5.,10.], pft_val=[6],
                    ensemble_path=ensemble_path,check_pixel=check_pixel )

        loc_lat   = [-26,-24]
        loc_lon   = [137.,138.]

        message   = hw_name +"_grass_wtd=10-15_"+time_spell

        plot_time_series_errorbar_select_regions( file_paths, var_name, date_s, date_e, seconds, hw_start,
                    loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name, message=message,
                    labal_names=labal_names, rain_val=10., wtd_val=[10.,15.], pft_val=[6],
                    ensemble_path=ensemble_path,check_pixel=check_pixel )


        loc_lat   = [-37,-36]
        loc_lon   = [147.,148.]
        message   = hw_name +"_forest_wtd=0-5_"+time_spell

        plot_time_series_errorbar_select_regions( file_paths, var_name, date_s, date_e, seconds, hw_start,
                    loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name, message=message,
                    labal_names=labal_names, rain_val=10., wtd_val=[0.,5.], pft_val=[2],
                    ensemble_path=ensemble_path,check_pixel=check_pixel )

        loc_lat   = [-31,-32]
        loc_lon   = [145.,146.]

        message   = hw_name +"_forest_wtd=5-10_"+time_spell

        plot_time_series_errorbar_select_regions( file_paths, var_name, date_s, date_e, seconds, hw_start,
                    loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name, message=message,
                    labal_names=labal_names, rain_val=10., wtd_val=[5.,10.], pft_val=[2],
                    ensemble_path=ensemble_path,check_pixel=check_pixel )

        loc_lat   = [-36,-34]
        loc_lon   = [141.,142.]

        message   = hw_name +"_forest_wtd=10-15_"+time_spell

        plot_time_series_errorbar_select_regions( file_paths, var_name, date_s, date_e, seconds, hw_start,
                    loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name, message=message,
                    labal_names=labal_names, rain_val=10., wtd_val=[10.,15.], pft_val=[2],
                    ensemble_path=ensemble_path,check_pixel=check_pixel )
