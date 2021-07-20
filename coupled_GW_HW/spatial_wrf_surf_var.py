#!/usr/bin/python

import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.animation as animation
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                        cartopy_ylim, latlon_coords, ALL_TIMES)

def plot_spatial_wrf_surf_press(file_path, case_name, timeidx):

    # Open the NetCDF file
    ncfile = Dataset(file_path)

    # Get the sea level pressure
    slp = getvar(ncfile, "slp",timeidx=timeidx)

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

    plt.title("Sea Level Pressure (hPa) of timestep-"+str(timeidx))

    fig.savefig("./plots/spatial_wrf_slp_"+case_name+"_"+str(timeidx), bbox_inches='tight', pad_inches=0.1)

    return ax

def plot_spatial_wrf_surf_var(file_path, case_name, var_name,val_min, var_max, timeidx):

    # Open the NetCDF file
    ncfile = Dataset(file_path)

    # Get the sea level pressure
    var = getvar(ncfile,var_name,timeidx=timeidx)

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

    plt.title("2m Temperature (K) of timestep-"+str(timeidx))

    fig.savefig("./plots/spatial_wrf_surf_"+var_name+"_"+case_name+"_"+str(timeidx), bbox_inches='tight', pad_inches=0.1)

    return ax

def plot_spatial_wrf_surf_var_diff(file_paths, case_names, var_name, val_min, var_max, timeidx):

    # for i, file_path in enumerate(file_paths):

    #     # Open the NetCDF file
    #     ncfile = Dataset(file_path)

    #     # Get the sea level pressure
    #     ('var_%d' % i) = getvar(ncfile,var_name,timeidx=timeidx)


    # Open the NetCDF file
    ncfile1 = Dataset(file_paths[0])
    ncfile2 = Dataset(file_paths[1])

    # Get the sea level pressure
    var_1 = getvar(ncfile1,var_name,timeidx=timeidx)
    var_2 = getvar(ncfile2,var_name,timeidx=timeidx)
    var_diff = var_2 - var_1
    # print(slp)

    # Smooth the sea level pressure since it tends to be noisy near the
    # mountains
    # smooth_var = smooth2d(var, 3, cenweight=4)

    # Get the latitude and longitude points
    lats, lons = latlon_coords(var_1)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(var_1)

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
    levels = np.arange(val_min, var_max+2, 2)
    # plt.contour(to_np(lons), to_np(lats), to_np(var_diff), levels=levels,  colors="black",
    #                     transform=crs.PlateCarree())
    plt.contourf(to_np(lons), to_np(lats), to_np(var_diff), levels=levels,
                         transform=crs.PlateCarree(), cmap=get_cmap("bwr")) #"jet" # "rainbow"

    # Add a color bar
    plt.colorbar(ax=ax, shrink=.98)

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(var_1))
    ax.set_ylim(cartopy_ylim(var_1))

    # Add the gridlines
    ax.gridlines(color="black", linestyle="dotted")

    plt.title("2m Temperature (K) of timestep-"+str(timeidx))

    fig.savefig("./plots/spatial_wrf_surf_diff_"+var_name+"_"+case_names[0]+"_vs_"+case_names[1]
                +"_"+str(timeidx), bbox_inches='tight', pad_inches=0.1)

    return ax

def plot_spatial_wrf_surf_var_diff_period_mean(file_paths, case_names, var_name,
                                        val_min, var_max, timeidx_s, timeidx_e, ts = None):

    # Open the NetCDF file
    ncfile1 = Dataset(file_paths[0])
    ncfile2 = Dataset(file_paths[1])

    # Get the variable
    var_tmp1 = getvar(ncfile1,var_name,timeidx=ALL_TIMES)
    var_tmp2 = getvar(ncfile2,var_name,timeidx=ALL_TIMES)
    # print(var_tmp1)
    if ts == None:
        var_1 = np.mean(var_tmp1[timeidx_s:timeidx_e],axis=0)
        var_2 = np.mean(var_tmp2[timeidx_s:timeidx_e],axis=0)
    else:
        print(var_tmp1[timeidx_s+ts:timeidx_e+ts:8])
        var_1 = np.mean(var_tmp1[timeidx_s+ts:timeidx_e+ts:8],axis=0)
        var_2 = np.mean(var_tmp2[timeidx_s+ts:timeidx_e+ts:8],axis=0)
    # print(var_1)
    var_diff = var_2 - var_1
    # print(slp)

    # Smooth the sea level pressure since it tends to be noisy near the
    # mountains
    # smooth_var = smooth2d(var, 3, cenweight=4)

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
    levels = np.arange(val_min, var_max+1, 1)
    # plt.contour(to_np(lons), to_np(lats), to_np(var_diff), levels=levels,  colors="black",
    #                     transform=crs.PlateCarree())
    plt.contourf(to_np(lons), to_np(lats), to_np(var_diff), levels=levels,
                         transform=crs.PlateCarree(), cmap=get_cmap("bwr")) #"jet" # "rainbow"

    # Add a color bar
    plt.colorbar(ax=ax, shrink=.98)

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(var_tmp1))
    ax.set_ylim(cartopy_ylim(var_tmp1))

    # Add the gridlines
    ax.gridlines(color="black", linestyle="dotted")

    # plt.title("2m Temperature (K) of ts="+str(ts))
    plt.title("2m Relative Humidity (%) of ts="+str(ts))
    if ts == None:
        fig.savefig("./plots/spatial_wrf_surf_diff_"+var_name+"_"+case_names[0]+"_vs_"+case_names[1]
                +"_period-"+str(timeidx_s)+"-"+str(timeidx_e), bbox_inches='tight', pad_inches=0.1)
    else:
        fig.savefig("./plots/spatial_wrf_surf_diff_"+var_name+"_"+case_names[0]+"_vs_"+case_names[1]
                +"_period-"+str(timeidx_s)+"-"+str(timeidx_e)+"_ts-"+str(ts), bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    print(sys.getrecursionlimit())
    sys.setrecursionlimit(50000)

    # # plot_spatial_wrf_surf_press
    # path      = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hires_r7264/WRF_output/"
    # file_name = "wrfout_d01_2012-12-01_00:00:00"
    # file_path = path+file_name
    # case_name = "GW"
    # ims = []
    # for timeidx in np.arange(0,249):
    #     ims.append([plot_spatial_wrf_surf_press(file_path,case_name,timeidx)])
    # print(ims)

    # # plot_spatial_wrf_surf_var
    # path      = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hires_r7264/WRF_output/"
    # file_name = "wrfout_d01_2012-12-01_00:00:00"
    # file_path = path+file_name
    # case_name = "GW"
    # var_name  = "T2"
    # val_min, var_max = 10, 50

    # ims = []
    # for timeidx in np.arange(0,249):
    #     ims.append([plot_spatial_wrf_surf_var(file_path, case_name, var_name, val_min, var_max, timeidx)])
    # print(ims)

    # ### plot_spatial_wrf_surf_var_diff
    # case_names = ['free_drain_hires_r7264','hires_r7264'] # the first case_name is set as control by default
    # file_name  = "wrfout_d01_2013-01-01_03:00:00"
    # var_name   = "T2"
    # val_min, var_max = -14, 14

    # file_paths = []
    # for case_name in case_names:
    #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/"
    #     file_path  = path + file_name
    #     file_paths.append(file_path)
    # print(file_paths)
    # ims = []
    # timeidxs = [3+8*3, 3+8*4, 3+8*5, 3+8*6, 3+8*7] # 12 pm at 4th - 8th Jan 2013
    # for timeidx in timeidxs:
    #     ims.append([plot_spatial_wrf_surf_var_diff(file_paths, case_names, var_name, val_min, var_max, timeidx)])

    ### plot_spatial_wrf_surf_var_diff_period_mean
    case_names = ['free_drain_11Jul','ctl_11Jul'] # the first case_name is set as control by default
    file_name  = "wrfout_d01_2013-01-01_03:00:00"
    var_name   = "rh2" #"T2"
    val_min, var_max = -10, 10

    file_paths = []
    for case_name in case_names:
        path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/"
        file_path  = path + file_name
        file_paths.append(file_path)

    timeidx_s = 0      # 12 pm at 4th - 8th Jan 2013
    timeidx_e = 14*8
    tss       = [ 0,   1,   2,    3,   4,   5,   6,    7 ]
              # 3am, 6am, 9am, 12pm, 3pm, 6pm, 9pm, 12am
              # 1pm, 4pm, 7pm, 10pm, 1am, 4am, 7am, 10am
    for ts in tss:
        plot_spatial_wrf_surf_var_diff_period_mean(file_paths, case_names, var_name, val_min, var_max, timeidx_s, timeidx_e, ts)

    # # Create the figure
    # fig = plt.figure(figsize=(12,6))
    # fig.subplots_adjust(hspace=0.3)
    # fig.subplots_adjust(wspace=0.2)

    # plt.rcParams['text.usetex']     = False
    # plt.rcParams['font.family']     = "sans-serif"
    # plt.rcParams['font.serif']      = "Helvetica"
    # plt.rcParams['axes.linewidth']  = 1.5
    # plt.rcParams['axes.labelsize']  = 14
    # plt.rcParams['font.size']       = 14
    # plt.rcParams['legend.fontsize'] = 12
    # plt.rcParams['xtick.labelsize'] = 12
    # plt.rcParams['ytick.labelsize'] = 14

    # almost_black = '#262626'
    # # change the tick colors also to the almost black
    # plt.rcParams['ytick.color']     = almost_black
    # plt.rcParams['xtick.color']     = almost_black

    # # change the text colors also to the almost black
    # plt.rcParams['text.color']      = almost_black

    # ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
    #       repeat_delay=1000)
    # ani.save("spatial_wrf_slp_"+case_name+".mp4")
