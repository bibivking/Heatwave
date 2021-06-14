#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature
from wrf import (getvar, interplevel, get_cartopy, cartopy_xlim,
                 cartopy_ylim, to_np, latlon_coords)

def plot_spatial_wrf_Tair_Wind(file_path,height,timeidx):

    ######################################################
    # Note that since I have updated Python, get_basemap cannot work 
    # for the new version. Thus the function cannot really work.
    ######################################################
    
    from wrf import get_basemap

    # Open the NetCDF file
    ncfile = Dataset(file_path)

    # Extract the pressure, geopotential height, and wind variables
    p = getvar(ncfile, "pressure",timeidx=timeidx)
    z = getvar(ncfile, "z", units="dm",timeidx=timeidx)
    ua = getvar(ncfile, "ua", units="kt",timeidx=timeidx)
    va = getvar(ncfile, "va", units="kt",timeidx=timeidx)
    wspd = getvar(ncfile, "wspd_wdir", units="kts",timeidx=timeidx)[0,:]

    # Interpolate geopotential height, u, and v winds to 500 hPa
    ht_hgt = interplevel(z, p, height)
    u_hgt = interplevel(ua, p, height)
    v_hgt = interplevel(va, p, height)
    wspd_hgt = interplevel(wspd, p, height)

    # Get the lat/lon coordinates
    lats, lons = latlon_coords(ht_hgt)

    # Get the basemap object
    bm = get_basemap(ht_hgt)

    # Create the figure
    fig = plt.figure(figsize=(12,9))
    ax = plt.axes()

    # Convert the lat/lon coordinates to x/y coordinates in the projection space
    x, y = bm(to_np(lons), to_np(lats))

    # Add the 500 hPa geopotential height contours
    levels = np.arange(520., 580., 6.)
    contours = bm.contour(x, y, to_np(ht_hgt), levels=levels, colors="black")
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i")

    # Add the wind speed contours
    levels = [25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 110, 120]
    wspd_contours = bm.contourf(x, y, to_np(wspd_hgt), levels=levels,
                                cmap=get_cmap("rainbow"))
    plt.colorbar(wspd_contours, ax=ax, orientation="horizontal", pad=.05)

    # Add the geographic boundaries
    bm.drawcoastlines(linewidth=0.25)
    bm.drawstates(linewidth=0.25)
    bm.drawcountries(linewidth=0.25)

    # Add the 500 hPa wind barbs, only plotting every 125th data point.
    bm.barbs(x[::125,::125], y[::125,::125], to_np(u_hgt[::125, ::125]),
            to_np(v_hgt[::125, ::125]), length=6)

    plt.title("500 MB Height (dm), Wind Speed (kt), Barbs (kt)")

    plt.show()

def plot_spatial_wrf(case_name,file_path,var_name,var_unit,height,timeidx,val_min,val_max):

    # Open the NetCDF file
    ncfile = Dataset(file_path)

    # Extract the pressure, geopotential height, and wind variables
    p    = getvar(ncfile, "pressure")
    var  = getvar(ncfile, var_name, units=var_unit,timeidx=timeidx)
    z    = getvar(ncfile, "z", units="dm",timeidx=timeidx)
    ua   = getvar(ncfile, "ua", units="m s-1",timeidx=timeidx)
    va   = getvar(ncfile, "va", units="m s-1",timeidx=timeidx)
    # wspd = getvar(ncfile, "wspd_wdir", units="m s-1",timeidx=timeidx)[0,:]

    # Interpolate geopotential height, u, and v winds to 500 hPa
    var_hgt  = interplevel(var, p, height)
    u_hgt    = interplevel(ua, p, height)
    v_hgt    = interplevel(va, p, height)
    z_hgt    = interplevel(z, p, 500)
    # wspd_hgt = interplevel(wspd, p, height)

    # Get the lat/lon coordinates
    lats, lons = latlon_coords(var_hgt)

    # Get the basemap object
    # bm         = get_basemap(var_hgt)
    # Get the cartopy mapping object
    cart_proj = get_cartopy(var_hgt)

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

    # Convert the lat/lon coordinates to x/y coordinates in the projection space
    # x, y = bm(to_np(lons), to_np(lats))

    # Add the 500 hPa geopotential height contours
    levels =  np.arange(512, 600, 4.)
    contours = plt.contour(to_np(lons), to_np(lats), to_np(z_hgt), levels = levels, colors="black",
                        transform=crs.PlateCarree())
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i")

    # Add the wind speed contours
    levels = np.arange(val_min, val_max, 5.)
    wspd_contours = plt.contourf(to_np(lons), to_np(lats), to_np(var_hgt), levels = levels,
                         transform=crs.PlateCarree(), cmap=get_cmap("coolwarm")) #"jet" #“rainbow”
    plt.colorbar(wspd_contours, ax=ax, orientation="horizontal", pad=.05)

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(var_hgt))
    ax.set_ylim(cartopy_ylim(var_hgt))

    # Add the 500 hPa wind barbs, only plotting every 125th data point.
    ax.quiver(to_np(lons[::5,::5]), to_np(lats[::5,::5]), to_np(u_hgt[::5, ::5]),
             to_np(v_hgt[::5, ::5]),  transform=crs.PlateCarree()) #scale=5,

    plt.title("500hPa Geopotential Height (dm), "+str(height)+"hPa Temperature (degC) and Barbs (m s-1) timestep-"+str(timeidx))

    fig.savefig("./plots/spatial_wrf_"+var_name+"_"+case_name+"_"+str(timeidx) , bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    path      = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/free_drain_hires_r7264/WRF_output/"
    file_name = "wrfout_d01_2012-12-01_00:00:00"
    file_path = path+file_name
    var_name  = "temp"
    var_unit  = "degC"
    height    = 850
    case_name = "FD"
    val_min, val_max = -50, 50
    for timeidx in np.arange(0,249):
      plot_spatial_wrf(case_name,file_path,var_name,var_unit,height,timeidx,val_min,val_max)
