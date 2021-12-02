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

def read_ensembles(file_path, var_name, time_s, time_e, loc_lat=None, loc_lon=None, mask_map=None):
    
    lat_name="XLAT"
    lon_name="XLONG"
    time = 
    if var_name == "Tair": 
        Tair = read_wrf_surf_var(file_path, "T2", loc_lat=loc_lat, loc_lon=loc_lon, mask_map=mask_map)
        tair = time_series_var(time,Tair,time_s,time_e)
    elif var_name == "PBL": 
        PBL  = read_wrf_surf_var(file_path, "PBLH", loc_lat=loc_lat, loc_lon=loc_lon, mask_map=mask_map)
        pbl  = time_series_var(time,PBL,time_s,time_e)
    elif var_name == "HC": 
        HC = calc_heat_content()
        hc  = time_series_var(time,HC,time_s,time_e)
    
    
    return tair, pbl, hc

def get_mask(land_file, time_s, time_e, rain_val=None, wtd_val=None, pft_val=None, check_pixel=False):
    
    # get pixels wants to keep in the whole domain --- Note that this function weren't tested
    
    if rain_val is not None:
        Time, rain= read_var(land_file, "Rainf_tavg", lat_name="lat", lon_name="lon")
        rain_mean = spital_var(UTC_to_AEST(Time), rain, time_s, time_e)
        rain_mask = np.ones(np.shape(rain_mean), dtype=bool)
        rain_mask = np.where( rain_mean*24.*3600. < rain_val, True, False)
        print("(rain_mask == True).sum()")
        print((rain_mask == True).sum())
        mask_map = rain_mask
        
    if wtd_val is not None:
        Time, wtd = read_var(land_file, "WaterTableD_tavg", lat_name="lat", lon_name="lon")
        wtd_mean = spital_var(UTC_to_AEST(Time), wtd, time_s, time_e)
        wtd_mask = np.ones(np.shape(wtd_mean), dtype=bool)
        wtd_mask = np.where(np.all([(wtd_mean/1000.) >=wtd_val[0], (wtd_mean/1000.) <wtd_val[1]],axis=0), True, False)
        print("(wtd_mask == True).sum()")
        print((wtd_mask == True).sum())
        if mask_map is not None:
            mask_map = np.all([mask_map,wtd_mask],axis=0)
        else:
            mask_map = wtd_mask
            
    if pft_val is not None:
        pft      = tree_mask(land_file,"Landcover_inst")
        pft_mask = np.ones(np.shape(pft), dtype=bool)
        if len(pft_val) == 1:
            pft_mask = np.where( pft == pft_val[0], True, False)
        elif (len(pft_val) == 2) :
            pft_mask = np.where( np.all([pft >=pft_val[0], pft<=pft_val[1]],axis=0), True, False))
        if mask_map is not None:
            mask_map = np.all([mask_map,pft_mask],axis=0)
        else:
            mask_map = pft_mask

    # check pixel cover
    if check_pixel:
    
        # check mask
        file_tmp = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg/wrfout_20090122-20090213_gw"
        file = nc.Dataset(file_tmp, mode='r')
        p1   = getvar(file, "pressure")
        lon  = file.variables['XLONG'][0,:,:]
        lat  = file.variables['XLAT'][0,:,:]
    
        fig1 = plt.figure(figsize=(12,9))
        # Get the lat/lon coordinates
    
        ax1 = plt.axes(projection=ccrs.PlateCarree())
        ax1.set_extent([110,155,-45,-10])
        ax1.coastlines(resolution="50m",linewidth=1)
    
        # Add gridlines
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top  = False
        gl.ylabels_right= False
        gl.xlines       = True
        gl.xlocator     = mticker.FixedLocator(np.arange(110,155,1))
        gl.ylocator     = mticker.FixedLocator(np.arange(-45,-10,1))
        gl.xformatter   = LONGITUDE_FORMATTER
        gl.yformatter   = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black','rotation': 90}
        gl.ylabel_style = {'size':10, 'color':'black'}
    
        mask = np.where(np.isnan(mask_map), 0, 1)
        plt.contourf(lon, lat, mask, levels=[-0.1,0.1,1.1], transform=ccrs.PlateCarree(), cmap=plt.cm.jet) 
        plt.savefig('./plots/figures/check_pixels_rain_wtd_pft.png',dpi=300)
        
    return mask_map

def plot_time_series( path, case_names, periods, time_ss, time_es, seconds=None,
                      loc_lats=None, loc_lons=None, message=None, rain_val=None, 
                      wtd_val=None, pft_val=None):
    
    print("======== In plot_time_series =========")

    # ============== read heatwave information ==============
    hw_num    = 0
    case_name = case_names[hw_num]
    period    = periods[hw_num]
    time_s    = time_ss[hw_num]
    time_e    = time_es[hw_num]
    loc_lat   = loc_lats[hw_num]
    loc_lon   = loc_lons[hw_num]
    
    if case_name == "hw2009_3Nov":
        rst_dates = ["20090117","20090118","20090119","20090120","20090121" ]
    elif case_name == "hw2013_3Nov":
        rst_dates = ["20121224","20121225","20121226","20121227","20121228" ]
    elif case_name == "hw2019_3Nov":
        rst_dates = ["20190103","20190104","20190105","20190106","20190107" ]
        
    # ============== set up file names ==============
    file_path  = path + case_name 
    file_paths = [ file_path + "/fd_rst_" + rst_dates[0] + "/WRF_output/wrfout_" + period,
                   file_path + "/fd_rst_" + rst_dates[1] + "/WRF_output/wrfout_" + period,
                   file_path + "/fd_rst_" + rst_dates[2] + "/WRF_output/wrfout_" + period,
                   file_path + "/fd_rst_" + rst_dates[3] + "/WRF_output/wrfout_" + period,
                   file_path + "/fd_rst_" + rst_dates[4] + "/WRF_output/wrfout_" + period,
                   file_path + "/gw_rst_" + rst_dates[0] + "/WRF_output/wrfout_" + period,
                   file_path + "/gw_rst_" + rst_dates[1] + "/WRF_output/wrfout_" + period,
                   file_path + "/gw_rst_" + rst_dates[2] + "/WRF_output/wrfout_" + period,
                   file_path + "/gw_rst_" + rst_dates[3] + "/WRF_output/wrfout_" + period,
                   file_path + "/gw_rst_" + rst_dates[4] + "/WRF_output/wrfout_" + period ]
    
    # ============== read variables ==============
    case_sum = len(file_paths)

    pbl      = np.zeros([case_sum, ntime])
    hc       = np.zeros([case_sum, ntime])
    tair     = np.zeros([case_sum, ntime])

    # ============== get mask =============
    is_mask     = (rain_val is not None) | (wtd_val is not None) | (pft_val is not None) 
    check_pixel = False 
    
    if is_mask:
        land_file  = path + case_name + "/ensemble_avg/LIS.CABLE."+period+"_gw.nc"
        print("land_file="+land_file)
        mask_map = get_mask(land_file, time_s, time_e, rain_val, wtd_val, pft_val, check_pixel=check_pixel)
        
    for case_num in np.arange(case_sum):
        file_path       = file_paths[case_num]
        tair[case_num,:]= read_ensembles(file_path, "Tair", time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, mask_map=mask_map)
        pbl[case_num,:] = read_ensembles(file_path, "PBL", time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, mask_map=mask_map)
        hc[case_num,:]  = read_ensembles(file_path, "HC", time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, mask_map=mask_map)

    # ============== plotting ==============
    fig, ax = plt.subplots()
    colors = cm.Paired(np.arange(0,12))

    delta = hw_start - datetime(2000,1,1)

    total_days = (time_e-time_s).days + 1

    values_mean = np.zeros((len(file_paths),total_days))
    values_max  = np.zeros((len(file_paths),total_days))
    values_min  = np.zeros((len(file_paths),total_days))



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

        var_mean, var_max, var_min = time_series_statistic(UTC_to_AEST(Time), Var_tmp, time_s, time_e) # shape = (selected_pixels_num, time_series)

        # ======== Local time 10am -4pm =======
        time = time_series_time(UTC_to_AEST(Time), time_s, time_e)

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
    # time2, var2 = time_series_var(Time2,Var2,time_s,time_e)
    # var = np.zeros((2,len(var1)))
    # var[0,:] = var1
    # var[1,:] = var2
    # ax.plot(t1, var*scale, alpha=0.5)
    # ax.set_ylabel('mm')
    ax.set_title(var_name)
    # ax.set_xticks(x1[::1440])
    # ax.set_xticklabels(np.arange(time_s,time_e,1))
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
    #        Set decks      #
    # #######################
    
    path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/" 
    
    case_names = ["hw2009_3Nov",
                  "hw2013_3Nov",
                  "hw2019_3Nov"]
    periods    = ["20090122-20090213",
                  "20121229-20130122",
                  "20190108-20190130"]
    time_ss    = [datetime(2009,1,28,0,0,0),
                  datetime(2013,1,4,0,0,0),
                  datetime(2019,1,14,0,0,0)]
    time_es    = [datetime(2009,2,13,23,59,0),
                  datetime(2013,1,22,23,59,0),
                  datetime(2019,1,30,23,59,0)]
    
    seconds       = [6.*60.*60.,18.*60.*60.]
    
    lat_name     = "lat"
    lon_name     = "lon"
    
    loc_lats   = [[-29,-28],
                  [-29,-28]]
                  
    loc_lons   = [[138.,139.],
                  [138.,139.]]
    
    
    if seconds= None:
        time_spell = "all"
    elif seconds[0] < seconds[1]:
        time_spell = "day"
    elif seconds[0] > seconds[1]:
        time_spell = "night"

    message   = time_spell

    plot_time_series(path, case_names, periods, time_ss, time_es, seconds,
                loc_lat=loc_lat, loc_lon=loc_lon, message=message,
                rain_val=10., wtd_val=[0.,5.], pft_val=[6],
                check_pixel=check_pixel )
