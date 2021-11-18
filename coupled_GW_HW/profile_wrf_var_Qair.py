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

def plot_profile_wrf(file_paths, var_name, var_units, time_s, time_e, seconds=None, message=None,calc_type=None):

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
        time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1,0,0,0))
    time = np.array(time_tmp)

    doy_tmp = [time[i].days for i in np.arange(len(time))]
    print(doy_tmp)

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

    doy = []
    for i in np.arange(len(time_cood)):
        doy.append(doy_tmp[i])


    # Get the WRF variables
    if var_name in ["temp","th","rh"]:
        z_tmp1      = getvar(ncfile1, "z", timeidx=ALL_TIMES)
        var_tmp1    = getvar(ncfile1, var_name, units=var_units, timeidx=ALL_TIMES)
    else:
        z_tmp1      = getvar(ncfile1, "z", timeidx=ALL_TIMES)
        print(z_tmp1)
        var_tmp1    = ncfile1.variables[var_name][:]
        print(var_tmp1)

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
    lat_slt = -36.
    lon_min = 139.0
    lon_max = 142.0

    start_point = CoordPair(lat=lat_slt, lon=lon_min) # (lat=-30., lon=115.0)
    end_point   = CoordPair(lat=lat_slt, lon=lon_max) # (lat=-30., lon=161.0)

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.
    var_out1    = []
    var_out2    = []

    for i in np.arange(np.shape(var1)[0]):
        var1_cross  = vertcross(var1[i], z1[i], wrfin=ncfile1, start_point=start_point,
                               end_point=end_point, latlon=True, meta=True)
        if i == 0:
            xy_loc     = np.linspace(lon_min,lon_max,len(var1_cross.coords['xy_loc']))
            vertical   = np.arange(0,5010.,10.)

            print(len(var1_cross.coords['xy_loc']))

        vertical_tmp   = to_np(var1_cross.coords['vertical'])[:]
        grid_x, grid_y = np.meshgrid(xy_loc,vertical_tmp)
        x              = np.reshape(grid_x,-1)
        y              = np.reshape(grid_y,-1)
        value          = np.reshape(to_np(var1_cross),-1)
        grid_X, grid_Y = np.meshgrid(xy_loc,vertical)

        var_out1.append(griddata((x, y), value, (grid_X, grid_Y), method="linear")[:])

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
            var_out2.append(griddata((x, y), value, (grid_X, grid_Y), method="linear")[:])

            grid_X, grid_Y, grid_x, grid_y = None, None, None, None
            x, y, value, var2_cross        = None, None, None, None

    print(var_out1)
    if calc_type == None:
        var_cross  = np.nanmean(var_out2,axis=0) - np.nanmean(var_out1,axis=0)
    elif calc_type == "Tmax":
        day_num = len(np.unique(doy))
        X_dim   = len(var_out1[0,:,0])
        Y_dim   = len(var_out1[0,0,:])
        tmax1   = np.zeros((X_dim,Y_dim))
        tmax2   = np.zeros((X_dim,Y_dim))
        for i in np.arange(day_num): 
            tmax1 = tmax1 + np.nanmax(var_out1[doy == doy[i],:,:],axis=0)
            tmax2 = tmax2 + np.nanmax(var_out2[doy == doy[i],:,:],axis=0)
        var_cross = (tmax2 - tmax1)/day_num
    elif calc_type == "Tmin":
        day_num = len(np.unique(doy))
        X_dim   = len(var_out1[0,:,0])
        Y_dim   = len(var_out1[0,0,:])
        tmin1   = np.zeros((X_dim,Y_dim))
        tmin2   = np.zeros((X_dim,Y_dim))
        for i in np.arange(day_num): 
            tmin1 = tmin1 + np.nanmin(var_out1[doy == doy[i],:,:],axis=0)
            tmin2 = tmin2 + np.nanmin(var_out2[doy == doy[i],:,:],axis=0)
        var_cross = (tmin2 - tmin1)/day_num


    # if len(file_paths) > 1:
    #     var_cross  = np.nanmean(var_out2,axis=0) - np.nanmean(var_out1,axis=0)
    # else:
    #     var_cross  = np.nanmean(var_out1,axis=0)

    print(np.shape(var_cross))

    # Create the figure that will have 3 subplots
    fig     = plt.figure(figsize=(12,9))
    ax      = fig.add_subplot(1,1,1)

    # Make the contour plot for var
    if var_name == "th":
        levels  = [-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.] 
    elif var_name == "rh":
        levels  = np.arange(-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,2,4,6,8,10,12,14,16,18,20)
    elif var_name == "temp":
        levels  = [-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.] 
    else:
        levels = None

    var_contours = ax.contourf( xy_loc, vertical, var_cross, levels=levels, cmap=get_cmap("coolwarm"))

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
    ax.set_xlabel("Longitude", fontsize=12)
    ax.set_ylabel("Geopotential Height (m)", fontsize=12)

    # Add titles

    if var_units == None:
        ax.set_title(var_name, {"fontsize" : 12})
    else:
        ax.set_title(var_name+" ("+var_units+")" , {"fontsize" : 12})

    if message == None:
        message = var_name
    else:
        message = message+"_"+var_name

    if calc_type == None:
        message = message
    else:
        message = message+"_"+calc_type

    fig.savefig("./plots/5Nov/wrf_prof/profile_wrf_"+message, bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":


    hw_name           = "hw2009_3Nov"
    var_name          = 'temp'
    var_unit          = 'degC'
    calc_type         = 'Tmax'
    # var_name          = 'th'
    # var_unit          = 'degC'

    seconds           = [6.*60.*60.,18.*60.*60.]

    # var_name          = 'QVAPOR' # "CLDFRA"
    # var_unit          = None

    # 2009
    if hw_name == "hw2009_3Nov":
        start_date= "20090122"
        end_date  = "20090213"
        time_s = datetime(2009,1,28,0,0,0,0)
        # time_e = datetime(2009,1,28,23,59,0,0)
        time_e = datetime(2009,2,8,23,59,0,0)
        # Time_s = datetime(2009,1,22,0,0,0,0)
        # Time_e = datetime(2009,2,13,23,59,0,0)

    elif hw_name == "hw2013_3Nov":
        start_date= "20121229"
        end_date  = "20130122"
        time_s = datetime(2013,1,4,0,0,0,0)
        time_e = datetime(2013,1,18,23,59,0,0)
        # Time_s = datetime(2012,12,29,0,0,0,0)
        # Time_e = datetime(2013,1,22,23,59,0,0)

    elif hw_name == "hw2019_3Nov":
        start_date= "20190108"
        end_date  = "20190130"
        time_s = datetime(2019,1,14,0,0,0)
        time_e = datetime(2019,1,26,23,59,0,0)
        # Time_s = datetime(2019,1,8,14,0,0,0)
        # Time_e = datetime(2019,1,30,0,0,0,0)
        #     
    cpl_atmo_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+hw_name+'/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_'+start_date+'-'+end_date+'_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_'+start_date+'-'+end_date+'_fd'  # atmo output of wrf-cable run

    file_paths        = [cpl_atmo_file_fd,cpl_atmo_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw

    if len(file_paths) > 1:
        message = "Daytime_GW-FD_"+str(time_s)+"-"+str(time_e)
    else:
        message = "Daytime_GW_"+str(time_s)+"-"+str(time_e)

    plot_profile_wrf(file_paths, var_name, var_unit, time_s, time_e, seconds=seconds, message=message,calc_type=calc_type)