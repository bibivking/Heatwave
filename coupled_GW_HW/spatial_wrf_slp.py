#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                         cartopy_ylim, latlon_coords)

def plot_spatial_wrf_slp(file_path,case_name,timeidx):

    # Open the NetCDF file
    ncfile = Dataset(file_path)

    # Get the sea level pressure
    slp = getvar(ncfile, "slp",timeidx=timeidx)

    print(slp)

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

    plt.title("Sea Level Pressure (hPa)")

    fig.savefig("./plots/spatial_wrf_slp_"+case_name+"_"+str(timeidx), bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    path      = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/free_drain_hires_r7264/WRF_output/"
    file_name = "wrfout_d01_2012-12-01_00:00:00"
    file_path = path+file_name
    case_name = "FD"
    for timeidx in np.arange(0,249):
        plot_spatial_wrf_slp(file_path,case_name,timeidx)
