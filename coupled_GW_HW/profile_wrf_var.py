#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from netCDF4 import Dataset
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (getvar, to_np, vertcross, smooth2d, CoordPair,
                 get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords)

def plot_profile_wrf_var(file_path, case_name, var_name, var_units, timeidx):

    # Open the NetCDF file
    ncfile     = Dataset(file_path)

    # Get the WRF variables
    slp        = getvar(ncfile, "slp", timeidx=timeidx)
    smooth_slp = smooth2d(slp, 3)
    z          = getvar(ncfile, "z", timeidx=timeidx)
    var        = getvar(ncfile, var_name, units=var_units, timeidx=timeidx)

    # Set the start point and end point for the cross section
    start_point = CoordPair(lat=-30., lon=115.0)
    end_point   = CoordPair(lat=-30., lon=160.0)

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.
    # z_cross    = vertcross(Z, z, wrfin=ncfile, start_point=start_point,
    #                     end_point=end_point, latlon=True, meta=True)
    var_cross = vertcross(var, z, wrfin=ncfile, start_point=start_point,
                           end_point=end_point, latlon=True, meta=True)
    print(var_cross)

    # Get the latitude and longitude points
    lats, lons = latlon_coords(slp)

    # Create the figure that will have 3 subplots
    fig     = plt.figure(figsize=(12,9))
    ax      = fig.add_subplot(1,1,1)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(slp)

    # Make the contour plot for var
    levels = np.arange(0,40,2)
    var_contours = ax.contourf(to_np(var_cross), levels=levels, cmap=get_cmap("jet"))
    print(var_contours)

    # Add the color bar
    cb_var = fig.colorbar(var_contours, ax=ax)
    cb_var.ax.tick_params(labelsize=5)

    # Set the x-ticks to use latitude and longitude labels.
    coord_pairs = to_np(var_cross.coords["xy_loc"])
    print(coord_pairs)

    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]
    print(x_labels)

    ax.set_xticks(x_ticks[::20])
    ax.set_xticklabels(x_labels[::20], rotation=45, fontsize=4)

    # Set the y-ticks to be height.
    vert_vals = to_np(var_cross.coords["vertical"])
    v_ticks = np.arange(vert_vals.shape[0])
    ax.set_yticks(v_ticks[::20])
    ax.set_yticklabels(vert_vals[::20], fontsize=4)

    # Set the x-axis and  y-axis labels
    ax.set_xlabel("Latitude, Longitude", fontsize=5)
    ax.set_ylabel("Height (m)", fontsize=5)

    # Add titles
    ax.set_title("Cross-Section of Wind Speed (kt)", {"fontsize" : 7})

    plt.show()

if __name__ == "__main__":

    ### plot_spatial_wrf_surf_var_diff_period_mean
    case_names = ['free_drain_hires_r7264','hires_r7264'] # the first case_name is set as control by default
    file_name  = "wrfout_d01_2013-01-01_03:00:00"
    var_name   = "th"
    var_units  = "degC"

    file_paths = []
    for case_name in case_names:
        path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/"
        file_path  = path + file_name
        file_paths.append(file_path)

    timeidx = 0      # 12 pm at 4th - 8th Jan 2013
    for file_path in file_paths:
        plot_profile_wrf_var(file_path, case_name, var_name, var_units, timeidx)
