#!/usr/bin/python

import sys
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.animation as animation
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                        cartopy_ylim, latlon_coords, ALL_TIMES)
from common_utils import *

def plot_spatial_wrf_surf_press(file_path, case_name, ts):

    # Open the NetCDF file
    ncfile = Dataset(file_path)

    # Get the sea level pressure
    slp = getvar(ncfile, "slp",timeidx=ts)

    # print(slp)

    # Smooth the sea level pressure since it tends to be noisy near the
    # mountains
    smooth_slp = smooth2d(slp, 3, cenweight=4)

    # Get the latitude and longitude points
    lats, lons = latlon_coords(slp)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(slp)

    # Create the figure
    fig = plt.figure(figsize=(12,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color']     = almost_black
    plt.rcParams['xtick.color']     = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']      = almost_black

    # Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=cart_proj)

    # Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                         facecolor="none",
                                         name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)

    # Make the contour outlines and filled contours for the smoothed sea level
    # pressure.
    plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp), 10, vmin=980, vmax=1020,  colors="black",
                        transform=crs.PlateCarree())
    plt.contourf(to_np(lons), to_np(lats), to_np(smooth_slp), 10, vmin=980, vmax=1020,
                         transform=crs.PlateCarree(), cmap=get_cmap("rainbow")) #"jet"

    # Add a color bar
    plt.colorbar(ax=ax, shrink=.98)

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(smooth_slp))
    ax.set_ylim(cartopy_ylim(smooth_slp))

    # Add the gridlines
    ax.gridlines(color="black", linestyle="dotted")

    plt.title("Sea Level Pressure (hPa) of timestep-"+str(ts))

    fig.savefig("./plots/spatial_wrf_slp_"+case_name+"_"+str(ts), bbox_inches='tight', pad_inches=0.1)

    return ax

def plot_spatial_wrf_surf_var_animation(file_path, case_name, var_name,val_min, var_max, ts):

    # Open the NetCDF file
    ncfile = Dataset(file_path)

    # Get the sea level pressure
    var = getvar(ncfile,var_name,timeidx=ts)

    # print(slp)

    # Smooth the sea level pressure since it tends to be noisy near the
    # mountains
    # smooth_var = smooth2d(var, 3, cenweight=4)

    # Get the latitude and longitude points
    lats, lons = latlon_coords(var)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(var)

    # Create the figure
    fig = plt.figure(figsize=(12,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color']     = almost_black
    plt.rcParams['xtick.color']     = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']      = almost_black

    # Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=cart_proj)

    # Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                         facecolor="none",
                                         name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)

    # Make the contour outlines and filled contours for the smoothed sea level
    # pressure.
    levels = np.arange(val_min, var_max, 2)
    plt.contour(to_np(lons), to_np(lats), to_np(var), levels=levels,  colors="black",
                        transform=crs.PlateCarree())
    plt.contourf(to_np(lons), to_np(lats), to_np(var), levels=levels,
                         transform=crs.PlateCarree(), cmap=get_cmap("rainbow")) #"jet"

    # Add a color bar
    plt.colorbar(ax=ax, shrink=.98)

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(var))
    ax.set_ylim(cartopy_ylim(var))

    # Add the gridlines
    ax.gridlines(color="black", linestyle="dotted")

    plt.title("2m Temperature (K) of timestep-"+str(ts))

    fig.savefig("./plots/spatial_wrf_surf_"+var_name+"_"+case_name+"_"+str(ts), bbox_inches='tight', pad_inches=0.1)

    return ax

def plot_spatial_wrf_surf_var(is_diff,file_paths, var_name, tss, message=None):

    # Open the NetCDF file
    ncfile1 = Dataset(file_paths[0])
    if is_diff:
        ncfile2 = Dataset(file_paths[1])

    for ts in tss:
        # Get the sea level pressure
        var1 = getvar(ncfile1, var_name, timeidx=ts)
        if is_diff:
            var2 = getvar(ncfile2, var_name, timeidx=ts)
            Var = var2 - var1
        else:
            Var = var1

        # Get the latitude and longitude points
        lats, lons = latlon_coords(var1)

        # Get the cartopy mapping object
        cart_proj = get_cartopy(var1)

        # Create the figure
        fig = plt.figure(figsize=(12,6))
        fig.subplots_adjust(hspace=0.3)
        fig.subplots_adjust(wspace=0.2)

        plt.rcParams['text.usetex']     = False
        plt.rcParams['font.family']     = "sans-serif"
        plt.rcParams['font.serif']      = "Helvetica"
        plt.rcParams['axes.linewidth']  = 1.5
        plt.rcParams['axes.labelsize']  = 14
        plt.rcParams['font.size']       = 14
        plt.rcParams['legend.fontsize'] = 12
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 14

        almost_black = '#262626'
        # change the tick colors also to the almost black
        plt.rcParams['ytick.color']     = almost_black
        plt.rcParams['xtick.color']     = almost_black

        # change the text colors also to the almost black
        plt.rcParams['text.color']      = almost_black

        # Set the GeoAxes to the projection used by WRF
        ax = plt.axes(projection=cart_proj)

        # Download and add the states and coastlines
        states = NaturalEarthFeature(category="cultural", scale="50m",
                                            facecolor="none",
                                            name="admin_1_states_provinces_shp")
        ax.add_feature(states, linewidth=.5, edgecolor="black")
        ax.coastlines('50m', linewidth=0.8)

        # Make the contour outlines and filled contours for the smoothed sea level
        # pressure.
        if is_diff:
            ranges = get_wrf_var_range_diff(var_name)
        else:
            ranges = None

        if get_reverse_colormap(var_name):
            cmap   = get_cmap("seismic_r")
        else:
            cmap   = get_cmap("seismic") # "bwr"

        if ranges == None:
            plt.contourf(to_np(lons), to_np(lats), to_np(Var),
                            transform=crs.PlateCarree(), cmap=cmap, extend='both')
        else:
            clevs = np.linspace(ranges[0],ranges[1], num=21)
            plt.contourf(to_np(lons), to_np(lats), to_np(Var), levels=clevs,
                            transform=crs.PlateCarree(), cmap=cmap, extend='both')

        # Add a color bar
        plt.colorbar(ax=ax, shrink=.98)

        # Set the map bounds
        ax.set_xlim(cartopy_xlim(var1))
        ax.set_ylim(cartopy_ylim(var1))

        # Add the gridlines
        ax.gridlines(color="black", linestyle="dotted")

        # time_name = ["10am", "1pm", "4pm", "7pm", "10pm", "1am", "4am", "7am"]
        # plt.title(var_name+" "+message+" day="+str(ts//8+1)+" "+time_name[ts%8])

        # if message == "201212":
            # time_name = ["10am", "1pm", "4pm", "7pm", "10pm", "1am", "4am", "7am"]
            # plt.title(var_name+" "+message+" day="+str(ts//8+1)+" "+time_name[ts%8])
        # elif message == "201301":
        #     time_name = ["1pm", "4pm", "7pm", "10pm", "1am", "4am", "7am", "10am"]
        #     plt.title(var_name+" "+message+" day="+str(ts//8+1)+" "+time_name[ts%8])
        if message == None:
            message = var_name
        else:
            message = message+"_"+var_name
        if is_diff:
            fig.savefig("./plots/5Nov/wrf_surf/spatial_wrf_surf_diff_"+message+"_"+str(ts), bbox_inches='tight', pad_inches=0.1)
        else:
            fig.savefig("./plots/5Nov/wrf_surf/spatial_wrf_surf_"+message+"_"+str(ts), bbox_inches='tight', pad_inches=0.1)

        Var = None

def plot_spatial_wrf_surf_var_period_mean(is_diff, file_paths, var_name, ts_s, ts_e, ts = None, message=None):

    # Open the NetCDF file
    ncfile1 = Dataset(file_paths[0])
    if is_diff:
        ncfile2 = Dataset(file_paths[1])

    # Get the variable
    var_tmp1 = getvar(ncfile1,var_name,timeidx=ALL_TIMES)
    if is_diff:
        var_tmp2 = getvar(ncfile2,var_name,timeidx=ALL_TIMES)

    if ts == None:
        var1 = np.mean(var_tmp1[ts_s:ts_e],axis=0)
        if is_diff:
            var2 = np.mean(var_tmp2[ts_s:ts_e],axis=0)
    else:
        print(var_tmp1[ts_s+ts:ts_e+ts:24])
        var1 = np.mean(var_tmp1[ts_s+ts:ts_e+ts:24],axis=0)
        if is_diff:
            var2 = np.mean(var_tmp2[ts_s+ts:ts_e+ts:24],axis=0)

    if is_diff:
        Var = var2 - var1
    else:
        Var = var1

    # Get the latitude and longitude points
    lats, lons = latlon_coords(var_tmp1)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(var_tmp1)

    # Create the figure
    fig = plt.figure(figsize=(12,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color']     = almost_black
    plt.rcParams['xtick.color']     = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']      = almost_black

    # Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=cart_proj)

    # Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                         facecolor="none",
                                         name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)

    # Make the contour outlines and filled contours for the smoothed sea level
    # pressure.
    if is_diff:
        ranges = get_wrf_var_range_diff(var_name)
    else:
        ranges = None

    if get_reverse_colormap(var_name):
        cmap   = get_cmap("bwr_r")
    else:
        cmap   = get_cmap("bwr")

    if ranges == None:
        plt.contourf(to_np(lons), to_np(lats), to_np(Var),
                        transform=crs.PlateCarree(), cmap=cmap, extend='both')
    else:
        clevs = np.linspace(ranges[0],ranges[1], num=21)
        plt.contourf(to_np(lons), to_np(lats), to_np(Var), levels=clevs,
                        transform=crs.PlateCarree(), cmap=cmap, extend='both')

    # Add a color bar
    plt.colorbar(ax=ax, shrink=.98)

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(var_tmp1))
    ax.set_ylim(cartopy_ylim(var_tmp1))

    # Add the gridlines
    ax.gridlines(color="black", linestyle="dotted")

    # plt.title("2m Temperature (K) of ts="+str(ts))
    # plt.title("2m Relative Humidity (%) of ts="+str(ts))
    plt.title(var_name+" ts="+str(ts))

    if message == None:
        message = var_name
    else:
        message = message + "_"+var_name

    if ts == None:
        fig.savefig("./plots/5Nov/wrf_surf/spatial_wrf_surf"+message
                +"_period-"+str(ts_s)+"-"+str(ts_e), bbox_inches='tight', pad_inches=0.1)
    else:
        fig.savefig("./plots/5Nov/wrf_surf/spatial_wrf_surf"+message
                +"_period-"+str(ts_s)+"-"+str(ts_e)+"_ts-"+str(ts), bbox_inches='tight', pad_inches=0.1)
    Var = None

def plot_spatial_wrf_surf(file_paths, var_name, time_s, time_e, loc_lat=None, loc_lon=None, message=None):

    # Open the NetCDF file
    encoding = 'utf-8' # Times in WRF output is btype, convert to string

    ncfile1  = Dataset(file_paths[0])
    ntime    = len(ncfile1.variables['Times'][:,0])
    time_tmp = []

    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(ncfile1.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(time_temp - datetime(2000,1,1))

    time = np.array(time_tmp)

    # to get lat and lon
    p1       = getvar(ncfile1, "pressure", timeidx=ALL_TIMES)

    # Extract the pressure, geopotential height, and wind variables
    Var1  = read_wrf_surf_var(file_paths[0], var_name, loc_lat, loc_lon)

    if var_name in ['T2','td2']:
        var1  = spital_var(time,Var1,time_s,time_e)-273.15
    else:
        scale = get_scale(var_name)
        var1  = spital_var(time,Var1,time_s,time_e)*scale


    if len(file_paths) > 1:
        Var2  = read_wrf_surf_var(file_paths[1], var_name, loc_lat, loc_lon)

        if var_name in ['T2','td2']:
            var2  = spital_var(time,Var2,time_s,time_e)-273.15
        else:
            scale = get_scale(var_name)
            var2  = spital_var(time,Var2,time_s,time_e)*scale

        # Calculate difference
        var = var2 - var1
    else:
        # Calculate difference
        var = var1

    # Get the lat/lon coordinates
    lats, lons = latlon_coords(p1)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(p1)

    # Create the figure
    fig = plt.figure(figsize=(12,9))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color']     = almost_black
    plt.rcParams['xtick.color']     = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']      = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor']  = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    # set the box type of sequence number
    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
    # choose colormap

    # Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=cart_proj)

    # Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                         facecolor="none",
                                         name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)

    # start plotting
    if loc_lat == None:
        ax.set_extent([135,155,-40,-25])
    else:
        ax.set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])


    # gaussian_filter(z,sigma=3)

    # Add the var contours
    if len(file_paths) > 1 and var_name == "T2":
        levels = np.linspace(-2., 2., num=20) # np.nanmin(var), np.nanmax(var)
    elif len(file_paths) > 1:
        max_val = np.abs(np.nanmax(var))
        min_val = np.abs(np.nanmin(var))
        max_range = np.maximum(max_val,min_val)
        levels = np.linspace(max_range*(-1.),max_range,num=20)
    else:
        levels = np.arange(np.nanmin(var), np.nanmax(var), 20)

        # levels = np.linspace(np.nanmin(var),np.nanmax(var), num=21) # np.nanmin(var), np.nanmax(var)
    var_contours = plt.contourf(to_np(lons), to_np(lats), to_np(var),
                   levels = levels, transform=crs.PlateCarree(), cmap=get_cmap("bwr"),extend='both') #,"jet" #“rainbow”#"coolwarm"
    plt.colorbar(var_contours, ax=ax, orientation="horizontal", pad=.05)

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(p1))
    ax.set_ylim(cartopy_ylim(p1))

    plt.title(var_name)

    if message == None:
        message = var_name
    else:
        message = message+"_"+var_name

    fig.savefig('./plots/5Nov/wrf_surf/15Oct/spatial_wrf_surf_'+message , bbox_inches='tight', pad_inches=0.1)


if __name__ == "__main__":

    #######################################################
    # Decks to run:
    #    plot_spital_map
    #######################################################

    var_3D = [
                'rh2',  # 2m Relative Humidity
                'T2',   # 2m Temperature
                'td2',  # 2m Dew Point Temperature
                'slp',  # Sea Level Pressure
                'ctt',  # Cloud Top Temperature
                'mdbz', # Maximum Reflectivity
                'pw',   # Precipitable Water
                'updraft_helicity', # Updraft Helicity
                'helicity',        # Storm Relative Helicity
                'cape_2d', # 2D CAPE (MCAPE/MCIN/LCL/LFC)
                'cloudfrac', # Cloud Fraction
              ]

                # 'ter',  # Model Terrain Height
    var_other         = ['RAINC','RAINNC','PSFC','U10','V10','TSK','PBLH']

    cpl_atmo_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_15Oct/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_20090122-20090213_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_20090122-20090213_fd'  # atmo output of wrf-cable run

    file_paths        = [cpl_atmo_file_fd,cpl_atmo_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw

    for j, var_name in enumerate(var_3D):

        # for i in np.arange(0,23):
        #     time_s = datetime(2009,1,22,0,0,0,0) + timedelta(days=int(i))
        #     time_e = datetime(2009,1,22,23,59,0,0) + timedelta(days=int(i))

        # 30 Jan
        for i in np.arange(0,22):
            time_s = datetime(2009,1,22,14,0,0,0) + timedelta(days=int(i))
            time_e = datetime(2009,1,23,13,59,0,0) + timedelta(days=int(i))

            if len(file_paths) > 1:
                message = "Couple_GW-FD_"+str(time_s)
            else:
                message = "Couple_GW_"+str(time_s)

            plot_spatial_wrf_surf(file_paths, var_name, time_s, time_e, message=message)


    #######################################################
    # Decks to run:
    #    plot_spatial_wrf_surf_press
    #    plot_spatial_wrf_surf_var
    #    plot_spatial_wrf_surf_var_period_mean
    #######################################################

    # if False:
    #     case_names = [  "hw2009_15Oct","hw2011_15Oct",
    #                     "hw2013_15Oct","hw2019_15Oct" ] # ["hw2014_15Oct","hw2017_15Oct"]

    #     file_names = [  "wrfout_20090122-20090213",
    #                     "wrfout_20110124-20110211",
    #                     "wrfout_20121226-20130114",
    #                     "wrfout_20190106-20190130" ]

    #     case_sum   = len(file_names)
    #     '''
    #     1. wrfout_d01_2012-12-01_00:00:00
    #     UTC   12am, 3am, 6am, 9am, 12pm, 3pm, 6pm, 9pm
    #     Local 10am, 1pm, 4pm, 7pm, 10pm, 1am, 4am, 7am
    #     2. wrfout_d01_2013-01-01_03:00:00
    #     UTC    3am, 6am, 9am, 12pm, 3pm, 6pm, 9pm, 12am
    #     Local  1pm, 4pm, 7pm, 10pm, 1am, 4am, 7am, 10am
    #     '''

    #     # ###################################
    #     #   plot_spatial_wrf_surf_var_diff  #
    #     # ###################################
    #     is_diff    = True
    #     var_name   = "rh2"

    #     for case_num in np.arange(case_sum):
    #         file_paths = []

    #         path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_names[case_num]+"/ensemble_avg/"
    #         file_path  = path + file_names[case_num]+"_fd"
    #         file_paths.append(file_path)
    #         file_path  = path + file_names[case_num]+"_gw"
    #         file_paths.append(file_path)

    #         tss = [0] #np.arange(0*8, 28*8, 8) # day 0-- day 28, local time 1pm

    #         message   = case_names[case_num]+"_GW-FD"
    #         plot_spatial_wrf_surf_var(is_diff, file_paths, var_name, tss, message=message)

    #     # ###############################################
    #     #   plot_spatial_wrf_surf_var_diff_period_mean  #
    #     # ###############################################
    #     var_name   = "rh2" #"rh2" #"T2"
    #     is_diff    = True #False

    #     ts_s       = [ 6*24, 6*24, 6*24, 6*24]
    #     ts_e       = [ 17*24, 13*24, 14*24, 20*24 ]
    #     for case_num in np.arange(case_sum):
    #         file_paths = []

    #         path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_names[case_num]+"/ensemble_avg/"
    #         file_path  = path + file_names[case_num]+"_fd"
    #         file_paths.append(file_path)
    #         file_path  = path + file_names[case_num]+"_gw"
    #         file_paths.append(file_path)

    #         message = case_names[case_num]+"_GW-FD"

    #         tss  = np.arange(24) # 1 pm
    #         for ts in tss:
    #             plot_spatial_wrf_surf_var_period_mean(is_diff,file_paths,var_name, ts_s[case_num], ts_e[case_num], ts, message)
