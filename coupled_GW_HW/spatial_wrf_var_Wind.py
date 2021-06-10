#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from wrf import getvar, interplevel, to_np, get_basemap, latlon_coords

def plot_spatial_wrf(case_name,file_path,var_name,var_unit,height,timeidx,val_min,val_max):

    # Open the NetCDF file
    ncfile = Dataset(file_path)

    # Extract the pressure, geopotential height, and wind variables
    p    = getvar(ncfile, "pressure")
    var  = getvar(ncfile, var_name, units=var_unit,timeidx=timeidx)
    # z    = getvar(ncfile, "z", units="dm")
    ua   = getvar(ncfile, "ua", units="m s-1")
    va   = getvar(ncfile, "va", units="m s-1")
    wspd = getvar(ncfile, "wspd_wdir", units="m s-1")[0,:]

    # Interpolate geopotential height, u, and v winds to 500 hPa
    var_hgt  = interplevel(var, p, height)
    u_hgt    = interplevel(ua, p, height)
    v_hgt    = interplevel(va, p, height)
    wspd_hgt = interplevel(wspd, p, height)

    # Get the lat/lon coordinates
    lats, lons = latlon_coords(var_hgt)

    # Get the basemap object
    bm         = get_basemap(var_hgt)

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

    ax = plt.axes()

    # Convert the lat/lon coordinates to x/y coordinates in the projection space
    x, y = bm(to_np(lons), to_np(lats))

    # Add the 500 hPa geopotential height contours
    levels = np.arange(val_min, val_max, 2.)
    contours = bm.contour(x, y, to_np(var_hgt), levels=levels, colors="black")
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i")

    # Add the wind speed contours
    levels = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    wspd_contours = bm.contourf(x, y, to_np(wspd_hgt), levels=levels,
                                cmap=get_cmap("rainbow"))
    plt.colorbar(wspd_contours, ax=ax, orientation="horizontal", pad=.05)

    # Add the geographic boundaries
    bm.drawcoastlines(linewidth=0.25)
    bm.drawstates(linewidth=0.25)
    bm.drawcountries(linewidth=0.25)

    # Add the 500 hPa wind barbs, only plotting every 125th data point.
    bm.barbs(x[::125,::125], y[::125,::125], to_np(u_hgt[::125, ::125]),
             to_np(v_hgt[::125, ::125]), length=5)

    plt.title(str(height)+" MB Height (dm), Wind Speed (kt), Barbs (kt)")

    fig.savefig("./plots/spatial_wrf_"+var_name+"_"+case_name+"_"+str(timeidx) , bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    path      = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/free_drain_hires_r7264/WRF_output/"
    file_name = "wrfout_d01_2012-12-01_00:00:00"
    file_path = path+file_name
    var_name  = "temp"
    var_unit  = "degC"
    height    = 850
    case_name = "FD"
    val_min, val_max = 0, 50
    for timeidx in np.arange(0,249):
      plot_spatial_wrf(case_name,file_path,var_name,var_unit,height,timeidx,val_min,val_max)
