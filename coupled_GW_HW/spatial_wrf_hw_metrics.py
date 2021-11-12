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

def plot_spatial_wrf_hw_metrics(file_paths, var_name, time_s, time_e, metric='Tmean',loc_lat=None, loc_lon=None, message=None):

    # Open the NetCDF file
    encoding = 'utf-8' # Times in WRF output is btype, convert to string

    ncfile1  = Dataset(file_paths[0])
    ntime    = len(ncfile1.variables['Times'][:,0])
    time_tmp = []

    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(ncfile1.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(time_temp - datetime(2000,1,1))

    time  = np.array(time_tmp)

    # to get lat and lon
    p1    = getvar(ncfile1, "pressure", timeidx=ALL_TIMES)

    # Extract the pressure, geopotential height, and wind variables
    Var1  = read_wrf_surf_var(file_paths[0], var_name, loc_lat, loc_lon)

    if metric == 'Tmean':
        var1  = spital_var(time,Var1,time_s,time_e)-273.15
    elif metric == 'Tmax':
        var1  = spital_var_max(time,Var1,time_s,time_e)-273.15

    if len(file_paths) > 1:
        Var2  = read_wrf_surf_var(file_paths[1], var_name, loc_lat, loc_lon)

        if metric == 'Tmean':
            var2  = spital_var(time,Var2,time_s,time_e)-273.15
        elif metric == 'Tmax':
            var2  = spital_var_max(time,Var2,time_s,time_e)-273.15

        var = var2 - var1
    else:
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
    if len(file_paths) > 1:
        # max_val = np.abs(np.nanmax(var))
        # min_val = np.abs(np.nanmin(var))
        # max_range = np.maximum(max_val,min_val)
        # levels = np.linspace(max_range*(-1.),max_range,num=20)

        if metric == 'Tmean':
            levels = np.arange(-4.0, 4.0, 0.4)
        elif metric == 'Tmax':
            levels = np.arange(-4.0, 4.0, 0.4)
    else:
        levels = np.arange(np.nanmin(var), np.nanmax(var), 21)

    var_contours = plt.contourf(to_np(lons), to_np(lats), to_np(var),
                   levels = levels[levels!=0], transform=crs.PlateCarree(), cmap=get_cmap("seismic"),extend='both') #,"jet" #“rainbow”#"coolwarm"
    plt.colorbar(var_contours, ax=ax, orientation="horizontal", pad=.05)  #"bwr"

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(p1))
    ax.set_ylim(cartopy_ylim(p1))

    plt.title(metric)

    if message == None:
        message = var_name+"_"+metric
    else:
        message = message+"_"+var_name+"_"+metric

    fig.savefig('./plots/5Nov/hw_metric/3Nov/spatial_wrf_hw_metrics_'+message , bbox_inches='tight', pad_inches=0.1)


if __name__ == "__main__":

    #######################################################
    # Decks to run:
    #    plot_spital_map
    #######################################################

    var_name = 'T2'

    cpl_atmo_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_20090122-20090213_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_20090122-20090213_fd'  # atmo output of wrf-cable run

    file_paths        = [cpl_atmo_file_fd,cpl_atmo_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw

    time_s = datetime(2009,1,28,14,0,0,0)
    time_e = datetime(2009,2,9,13,59,0,0)

    if len(file_paths) > 1:
        message = 'Couple_GW-FD_'+str(time_s)+'-'+str(time_e)
    else:
        message = 'Couple_GW_'+str(time_s)+'-'+str(time_e)

    metric   = 'Tmean'
    plot_spatial_wrf_hw_metrics(file_paths, var_name, time_s, time_e, metric=metric, message=message)

    metric   = 'Tmax'
    plot_spatial_wrf_hw_metrics(file_paths, var_name, time_s, time_e, metric=metric, message=message)
