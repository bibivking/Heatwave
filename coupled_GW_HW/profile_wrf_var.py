#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from netCDF4 import Dataset
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (getvar, to_np, vertcross, CoordPair,
                 get_cartopy, latlon_coords, ALL_TIMES)

def plot_profile_wrf_var(file_path, case_name, var_name, var_units, timeidx):

    # Open the NetCDF file
    ncfile     = Dataset(file_path)

    # Get the WRF variables
    z          = getvar(ncfile, "z", timeidx=timeidx)
    print(z)
    var        = getvar(ncfile, var_name, units=var_units, timeidx=timeidx)

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

    # Add titles
    ax.set_title(var_name+" ("+var_units+") ts-"+str(timeidx), {"fontsize" : 12})

    fig.savefig("./plots/profile_wrf_"+var_name+"_"+case_name+"_"+str(timeidx) , bbox_inches='tight', pad_inches=0.1)

def plot_profile_wrf_var_diff_period_mean(file_paths, case_names, var_name, var_units, timeidx_s, timeidx_e, ts=None):

    # Open the NetCDF file
    ncfile1     = Dataset(file_paths[0])
    ncfile2     = Dataset(file_paths[1])

    # Get the WRF variables
    z           = getvar(ncfile1, "z")
    var_tmp1    = getvar(ncfile1, var_name, units=var_units, timeidx=ALL_TIMES)
    var_tmp2    = getvar(ncfile2, var_name, units=var_units, timeidx=ALL_TIMES)    

    if ts == None:
        var_1 = np.mean(var_tmp1[timeidx_s:timeidx_e],axis=0)
        var_2 = np.mean(var_tmp2[timeidx_s:timeidx_e],axis=0)
    else:
        print(var_tmp1[timeidx_s+ts:timeidx_e+ts:8])
        var_1 = np.mean(var_tmp1[timeidx_s+ts:timeidx_e+ts:8],axis=0)
        var_2 = np.mean(var_tmp2[timeidx_s+ts:timeidx_e+ts:8],axis=0)

    var         = var_2 - var_1

    # Set the start point and end point for the cross section
    start_point = CoordPair(lat=-30., lon=115.0)
    end_point   = CoordPair(lat=-30., lon=161.0)

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.
    var_cross  = vertcross(var,  z, wrfin=ncfile1, start_point=start_point,
                           end_point=end_point, latlon=True, meta=True)

    # Create the figure that will have 3 subplots
    fig     = plt.figure(figsize=(12,9))
    ax      = fig.add_subplot(1,1,1)

    # Make the contour plot for var
    levels  = np.arange(-5, 6, 1)
    var_contours = ax.contourf(to_np(var_cross), levels=levels, cmap=get_cmap("coolwarm"))

    # Add the color bar
    cb_var = fig.colorbar(var_contours, ax=ax)
    cb_var.ax.tick_params(labelsize=12)

    # Set the x-ticks to use latitude and longitude labels.
    coord_pairs = to_np(var_cross.coords["xy_loc"])

    x_ticks     = np.arange(coord_pairs.shape[0])
    x_labels    = [pair.latlon_str() for pair in to_np(coord_pairs)]

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

    # Add titles
    ax.set_title(var_name+" ("+var_units+") ts-"+str(ts), {"fontsize" : 12})

    if ts == None:
        fig.savefig("./plots/profile_wrf_diff_"+var_name+"_"+case_names[0]+"_vs_"+case_names[1]
                +"_period-"+str(timeidx_s)+"-"+str(timeidx_e), bbox_inches='tight', pad_inches=0.1)
    else:
        fig.savefig("./plots/profile_wrf_diff_"+var_name+"_"+case_names[0]+"_vs_"+case_names[1]
                +"_period-"+str(timeidx_s)+"-"+str(timeidx_e)+"_ts-"+str(ts), bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    ### plot_profile_wrf_var
    case_names = ['free_drain_hires_r7264','hires_r7264'] # the first case_name is set as control by default
    file_name  = "wrfout_d01_2013-01-01_03:00:00"
    var_name   = "th"
    var_units  = "degC"

    file_paths = []
    for case_name in case_names:
        path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/"
        file_path  = path + file_name
        file_paths.append(file_path)


    timeidx_s = 0      # 12 pm at 4th - 8th Jan 2013
    timeidx_e = 14*8
    ts        = 3      # 12 pm
    plot_profile_wrf_var_diff_period_mean(file_paths, case_names, var_name, var_units, timeidx_s, timeidx_e, ts)