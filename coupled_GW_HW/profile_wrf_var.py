#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from netCDF4 import Dataset
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from scipy.interpolate import griddata
from wrf import (getvar, to_np, vertcross, CoordPair,
                 get_cartopy, latlon_coords, ALL_TIMES)
from common_utils import *

def plot_profile_wrf_var(file_path, case_name, var_name, var_units, ts, message = None):

    # Open the NetCDF file
    ncfile     = Dataset(file_path)

    # Get the WRF variables
    z          = getvar(ncfile, "z", timeidx=ts)
    print(z)
    if var_units == None:
        var    = getvar(ncfile, var_name, timeidx=ts)
    else:
        var    = getvar(ncfile, var_name, units=var_units, timeidx=ts)

    # Set the start point and end point for the cross section
    start_point = CoordPair(lat=-30., lon=115.0)
    end_point   = CoordPair(lat=-30., lon=161.0)

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.
    var_cross = vertcross(var, z, wrfin=ncfile, start_point=start_point,
                           end_point=end_point, latlon=True, meta=True)

    # Get the latitude and longitude points
    lats, lons = latlon_coords(var)

    # Create the figure that will have 3 subplots
    fig     = plt.figure(figsize=(12,9))
    ax      = fig.add_subplot(1,1,1)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(var)

    # Make the contour plot for var
    levels = np.arange(0, 50,2)
    var_contours = ax.contourf(to_np(var_cross), levels=levels, cmap=get_cmap("jet"))

    # Add the color bar
    cb_var = fig.colorbar(var_contours, ax=ax)
    cb_var.ax.tick_params(labelsize=12)

    # Set the x-ticks to use latitude and longitude labels.
    coord_pairs = to_np(var_cross.coords["xy_loc"])

    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]

    ax.set_xticks(x_ticks[::20])
    ax.set_xticklabels(x_labels[::20], rotation=45, fontsize=12)

    # Set the y-ticks to be height.
    vert_vals = to_np(var_cross.coords["vertical"])
    v_ticks = np.arange(vert_vals.shape[0])
    ax.set_yticks(v_ticks[::20])
    ax.set_yticklabels(vert_vals[::20], fontsize=12)

    # Set the x-axis and  y-axis labels
    ax.set_xlabel("Latitude, Longitude", fontsize=12)
    ax.set_ylabel("Height (m)", fontsize=12)

    if message == None:
        message = var_name
    else:
        message = message + "_" + var_name
    # Add titles
    if var_units == None:
        ax.set_title(var_name+" ts-"+str(ts), {"fontsize" : 12})
    else:
        ax.set_title(var_name+" ("+var_units+") ts-"+str(ts), {"fontsize" : 12})
    fig.savefig("./plots/5Nov/wrf_prof/profile_wrf_"+message+"_ts-"+str(ts) , bbox_inches='tight', pad_inches=0.1)

def plot_profile_wrf_var_diff_period_mean(file_paths, var_name, var_units, ts_s, ts_e, ts=None, message=None):

    # Open the NetCDF file
    ncfile1     = Dataset(file_paths[0])
    ncfile2     = Dataset(file_paths[1])

    # Get the WRF variables
    z_temp1     = getvar(ncfile1, "z", timeidx=ALL_TIMES)
    z_temp2     = getvar(ncfile2, "z", timeidx=ALL_TIMES)

    if var_units == None:
        var_tmp1    = getvar(ncfile1, var_name, timeidx=ALL_TIMES)
        var_tmp2    = getvar(ncfile2, var_name, timeidx=ALL_TIMES)
    else:
        var_tmp1    = getvar(ncfile1, var_name, units=var_units, timeidx=ALL_TIMES)
        var_tmp2    = getvar(ncfile2, var_name, units=var_units, timeidx=ALL_TIMES)

    if ts == None:
        var_1 = np.mean(var_tmp1[ts_s:ts_e],axis=0)
        var_2 = np.mean(var_tmp2[ts_s:ts_e],axis=0)
        z_1   = np.mean(z_temp1[ts_s:ts_e],axis=0)
        z_2   = np.mean(z_temp2[ts_s:ts_e],axis=0)
    else:
        var_1 = np.mean(var_tmp1[ts_s+ts:ts_e:24],axis=0)
        var_2 = np.mean(var_tmp2[ts_s+ts:ts_e:24],axis=0)
        z_1   = np.mean(z_temp1[ts_s+ts:ts_e:24],axis=0)
        z_2   = np.mean(z_temp2[ts_s+ts:ts_e:24],axis=0)

    # Set the start point and end point for the cross section
    # start_point = CoordPair(lat=-30., lon=115.0)
    # end_point   = CoordPair(lat=-30., lon=161.0)
    start_point = CoordPair(lat=-30., lon=135.0)
    end_point   = CoordPair(lat=-30., lon=155.0)

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.
    var1_cross  = vertcross(var_1, z_1, wrfin=ncfile1, start_point=start_point,
                           end_point=end_point, latlon=True, meta=True)
    var2_cross  = vertcross(var_2, z_2, wrfin=ncfile2, start_point=start_point,
                           end_point=end_point, latlon=True, meta=True)

    print(var1_cross.coords['vertical'][0:10])
    print(var2_cross.coords['vertical'][0:10])
    var_cross  = to_np(var2_cross) - to_np(var1_cross)

    # Create the figure that will have 3 subplots
    fig     = plt.figure(figsize=(12,9))
    ax      = fig.add_subplot(1,1,1)

    # Make the contour plot for var
    if var_name == "th":
        levels  = np.arange(-1, 1, 0.2)
    elif var_name == "rh":
        levels  = np.arange(-20, 20, 1)
    elif var_name == "temp":
        levels  = np.arange(-1., 1., 0.2)
    print(np.shape(var_cross))
    var_contours = ax.contourf( var1_cross.coords["vertical"], var1_cross.coords["xy_loc"], var_cross, levels=levels, cmap=get_cmap("coolwarm"))

    # Add the color bar
    cb_var = fig.colorbar(var_contours, ax=ax)
    cb_var.ax.tick_params(labelsize=12)

    # # Set the x-ticks to use latitude and longitude labels.
    # coord_pairs = to_np(var1_cross.coords["xy_loc"])
    # x_ticks     = np.arange(coord_pairs.shape[0])
    # x_labels    = [pair.latlon_str() for pair in to_np(coord_pairs)]
    # ax.set_xticks(x_ticks[::20])
    # ax.set_xticklabels(x_labels[::20], rotation=45, fontsize=12)

    # # Set the y-ticks to be height.
    # vert_vals = to_np(var1_cross.coords["vertical"])
    # v_ticks   = np.arange(vert_vals.shape[0])
    # ax.set_yticks(vert_vals)
    # ax.set_yticklabels(vert_vals, fontsize=12)
    # ax.set_ylim((0,2000))

    # Set the x-axis and  y-axis labels
    ax.set_xlabel("Latitude, Longitude", fontsize=12)
    ax.set_ylabel("Geopotential Height (m)", fontsize=12)

    # Add titles

    if var_units == None:
        ax.set_title(var_name+" ts-"+str(ts), {"fontsize" : 12})
    else:
        ax.set_title(var_name+" ("+var_units+") ts-"+str(ts), {"fontsize" : 12})

    if message == None:
        message = var_name
    else:
        message = message+"_"+var_name

    if ts == None:
        fig.savefig("./plots/5Nov/wrf_prof/profile_wrf_"+message
                +"_tss-"+str(ts_s)+"-"+str(ts_e), bbox_inches='tight', pad_inches=0.1)
    else:
        fig.savefig("./plots/5Nov/wrf_prof/profile_wrf_"+message
                +"_tss-"+str(ts_s)+"-"+str(ts_e)+"_ts-"+str(ts), bbox_inches='tight', pad_inches=0.1)

def plot_profile_wrf(file_paths, var_name, var_units, time_s, time_e, seconds=None, message=None):

    # Open the NetCDF file
    ncfile1     = Dataset(file_paths[0])

    # process time
    encoding    = 'utf-8' # Times in WRF output is btype, convert to string
    Time_s      = time_s - datetime(2000,1,1,0,0,0)
    Time_e      = time_e - datetime(2000,1,1,0,0,0)
    ntime       = len(ncfile1.variables['Times'][:,0])
    time_tmp    = []

    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(ncfile1.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(time_temp - datetime(2000,1,1,0,0,0))
    time = np.array(time_tmp)

    if seconds == None:
        time_cood = (time>=Time_s) & (time<Time_e)
    else:
        time_cood = []
        for j in np.arange(len(time)):
            if seconds[0] >= seconds[1]:
                if_seconds = (time[j].seconds >= seconds[0]) | (time[j].seconds < seconds[1])
            else:
                if_seconds = (time[j].seconds >= seconds[0]) & (time[j].seconds < seconds[1])
            time_cood.append((time[j]>=Time_s) & (time[j]<Time_e) & if_seconds)

    # Get the WRF variables
    z_tmp1      = getvar(ncfile1, "z", timeidx=ALL_TIMES)
    var_tmp1    = getvar(ncfile1, var_name, units=var_units, timeidx=ALL_TIMES)
    z1          = z_tmp1[time_cood,:,:,:]
    var1        = var_tmp1[time_cood,:,:,:]

    if len(file_paths) > 1:
        ncfile2     = Dataset(file_paths[1])
        z_tmp2      = getvar(ncfile2, "z", timeidx=ALL_TIMES)
        var_tmp2    = getvar(ncfile2, var_name, units=var_units, timeidx=ALL_TIMES)
        z2          = z_tmp2[time_cood,:,:,:]
        var2        = var_tmp2[time_cood,:,:,:]
        print(np.shape(var2))

    # Set the start point and end point for the cross section
    start_point = CoordPair(lat=-36., lon=135.0) # (lat=-30., lon=115.0)
    end_point   = CoordPair(lat=-36., lon=155.0) # (lat=-30., lon=161.0)

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.
    var_out1    = []
    var_out2    = []

    for i in np.arange(np.shape(var1)[0]):
        var1_cross  = vertcross(var1[i], z1[i], wrfin=ncfile1, start_point=start_point,
                               end_point=end_point, latlon=True, meta=True)
        if i == 0:
            xy_loc     = var1_cross.coords['xy_loc'][:]
            vertical   = np.arange(0,5010.,10.)

            print("======var1_cross.coords['vertical']=====")
            print(xy_loc)

        vertical_tmp   = to_np(var1_cross.coords['vertical'])[:]
        grid_x, grid_y = np.meshgrid(xy_loc,vertical_tmp)
        x              = np.reshape(grid_x,-1)
        y              = np.reshape(grid_y,-1)
        value          = np.reshape(to_np(var1_cross),-1)
        grid_X, grid_Y = np.meshgrid(xy_loc,vertical)
        #
        # print("======grid_x=====")
        # print(grid_x)
        #
        # print("======x=====")
        # print(x)
        #
        # print("======value=====")
        # print(value)
        #
        # print("======grid_X=====")
        # print(grid_X)
        var_out1.append(griddata((x, y), value, (grid_X, grid_Y), method="linear"))

        grid_X, grid_Y, grid_x, grid_y        = None, None, None, None
        x, y, value, vertical_tmp, var1_cross = None, None, None, None, None

        if len(file_paths) > 1:
            var2_cross  = vertcross(var2[i], z2[i], wrfin=ncfile2, start_point=start_point,
                                   end_point=end_point, latlon=True, meta=True)
            vertical_tmp   = to_np(var2_cross.coords['vertical'])[:]
            grid_x, grid_y = np.meshgrid(xy_loc,vertical_tmp)
            x              = np.reshape(grid_x,-1)
            y              = np.reshape(grid_y,-1)
            value          = np.reshape(to_np(var2_cross),-1)
            grid_X, grid_Y = np.meshgrid(xy_loc,vertical)
            var_out2.append(griddata((x, y), value, (grid_X, grid_Y), method="linear"))

            grid_X, grid_Y, grid_x, grid_y = None, None, None, None
            x, y, value, var2_cross        = None, None, None, None

    if len(file_paths) > 1:
        var_cross  = np.nanmean(var2_cross,axis=0) - np.nanmean(var1_cross,axis=0)
    else:
        var_cross  = np.nanmean(var1_cross,axis=0)

    print(np.shape(var_cross))

    # Create the figure that will have 3 subplots
    fig     = plt.figure(figsize=(12,9))
    ax      = fig.add_subplot(1,1,1)

    # Make the contour plot for var
    if var_name == "th":
        levels  = np.arange(-1, 1, 0.2)
    elif var_name == "rh":
        levels  = np.arange(-20, 20, 1)
    elif var_name == "temp":
        levels  = np.arange(-1., 1., 0.2)
    var_contours = ax.contourf( vertical, xy_loc, var_cross, levels=levels, cmap=get_cmap("coolwarm"))

    # Add the color bar
    cb_var = fig.colorbar(var_contours, ax=ax)
    cb_var.ax.tick_params(labelsize=12)

    # # Set the x-ticks to use latitude and longitude labels.
    # coord_pairs = to_np(var1_cross.coords["xy_loc"])
    # x_ticks     = np.arange(coord_pairs.shape[0])
    # x_labels    = [pair.latlon_str() for pair in to_np(coord_pairs)]
    # ax.set_xticks(x_ticks[::20])
    # ax.set_xticklabels(x_labels[::20], rotation=45, fontsize=12)

    # # Set the y-ticks to be height.
    # vert_vals = to_np(var1_cross.coords["vertical"])
    # v_ticks   = np.arange(vert_vals.shape[0])
    # ax.set_yticks(vert_vals)
    # ax.set_yticklabels(vert_vals, fontsize=12)
    # ax.set_ylim((0,2000))

    # Set the x-axis and  y-axis labels
    ax.set_xlabel("Latitude, Longitude", fontsize=12)
    ax.set_ylabel("Geopotential Height (m)", fontsize=12)

    # Add titles

    if var_units == None:
        ax.set_title(var_name+" ts-"+str(ts), {"fontsize" : 12})
    else:
        ax.set_title(var_name+" ("+var_units+") ts-"+str(ts), {"fontsize" : 12})

    if message == None:
        message = var_name
    else:
        message = message+"_"+var_name

    fig.savefig("./plots/5Nov/wrf_prof/profile_wrf_"+message, bbox_inches='tight', pad_inches=0.1)


if __name__ == "__main__":

    cpl_atmo_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_20090122-20090213_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_20090122-20090213_fd'  # atmo output of wrf-cable run

    file_paths        = [cpl_atmo_file_fd,cpl_atmo_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw

    var_name          = 'th'
    var_unit          = 'degC'

    seconds           = [20.*60.*60.,8.*60.*60.]
    time_s            = datetime(2009,1,29,14,0,0,0)
    time_e            = datetime(2009,1,30,13,59,0,0)

    if len(file_paths) > 1:
        message = "Daytime_GW-FD_"+str(time_s)
    else:
        message = "Daytime_GW_"+str(time_s)

    plot_profile_wrf(file_paths, var_name, var_unit, time_s, time_e, seconds=seconds, message=message)

    # ### ========================= Old Deck ============================
    # case_names = [  "hw2009_15Oct", "hw2011_15Oct",
    #                 "hw2013_15Oct", "hw2019_15Oct" ] # ["hw2014_15Oct","hw2017_15Oct"]
    #
    # file_names = [  "wrfout_20090122-20090213",
    #                 "wrfout_20110124-20110211",
    #                 "wrfout_20121226-20130114",
    #                 "wrfout_20190106-20190130" ]
    #
    # case_sum   = len(file_names)
    #
    # ts_s       = [ 6*24, 6*24, 6*24, 6*24]
    # ts_e       = [ 17*24, 13*24, 14*24, 20*24 ]
    #
    # ##########################################
    # # plot_profile_wrf_var_diff_period_mean  #
    # ##########################################
    #
    # var_name   = 'temp' #"th" #"rh"
    # var_units  = "degC" #None #'%' #"degC"
    #
    # # for case_num in np.arange(case_sum):
    # case_num   = 0
    # file_paths = []
    #
    # path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_names[case_num]+"/ensemble_avg/"
    # file_path  = path + file_names[case_num]+"_fd"
    # file_paths.append(file_path)
    # file_path  = path + file_names[case_num]+"_gw"
    # file_paths.append(file_path)
    #
    # tss       = np.arange(24) # 1 pm
    # message   = case_names[case_num]+"_GW-FD"
    # # for ts in tss:
    # ts = 0
    # plot_profile_wrf_var_diff_period_mean(file_paths, var_name, var_units, ts_s[case_num], ts_e[case_num], ts, message=message)
    #
    # # var_name   = 'rh' #"th" #"rh"
    # # var_units  = None #None #'%' #"degC"
    #
    # # for case_num in np.arange(case_sum):
    # #     file_paths = []
    #
    # #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_names[case_num]+"/ensemble_avg/"
    # #     file_path  = path + file_names[case_num]+"_fd"
    # #     file_paths.append(file_path)
    # #     file_path  = path + file_names[case_num]+"_gw"
    # #     file_paths.append(file_path)
    #
    # #     tss       = np.arange(24) # 1 pm
    # #     message   = case_names[case_num]+"_GW-FD"
    # #     for ts in tss:
    # #         plot_profile_wrf_var_diff_period_mean(file_paths, var_name, var_units, ts_s[case_num], ts_e[case_num], ts, message=message)
