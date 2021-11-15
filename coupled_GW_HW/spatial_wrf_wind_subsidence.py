#!/usr/bin/python

import sys
import numpy as np
import pint
from netCDF4 import Dataset,num2date
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature
import metpy.calc as mpcalc
from metpy.units import units
import scipy.ndimage as ndimage
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                        cartopy_ylim, latlon_coords, ALL_TIMES)
from common_utils import *

def wind_subsidence(file_paths, var_name, var_unit, height, time_s, time_e, loc_lat=None, loc_lon=None, message=None):

    # Open the NetCDF file
    encoding = 'utf-8' # Times in WRF output is btype, convert to string

    ncfile1  = Dataset(file_paths[0])
    ntime    = len(ncfile1.variables['Times'][:,0])
    lat      = ncfile1.variables['XLAT'][0,:,:]
    lon      = ncfile1.variables['XLONG'][0,:,:]

    time_tmp = []

    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(ncfile1.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1))

    time = np.array(time_tmp)

    # to get lat and lon
    p1    = getvar(ncfile1, "pressure", timeidx=ALL_TIMES)

    Var1  = read_wrf_hgt_var(file_paths[0], var_name, var_unit, height, loc_lat, loc_lon)# "wa", "m s-1",
    var1  = spital_var(time,Var1,time_s,time_e)

    if len(file_paths) > 1:
        Var2  = read_wrf_hgt_var(file_paths[1], var_name, var_unit, height, loc_lat, loc_lon) # "wa", "m s-1",
        var2  = spital_var(time,Var2,time_s,time_e)
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

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(p1))
    ax.set_ylim(cartopy_ylim(p1))

    # ------- Plot wa -------
    if len(file_paths) > 1:
        levels = [-30., -25., -20., -15., -10., -5., 5., 10., 15., 20., 25., 30.]
        # np.arange(-0.02, 0.021, 0.002)
    else:
        levels = np.arange(-100., 100., 5.)

    if var_name =="omg":
        cf = ax.contourf(to_np(lons), to_np(lats), var, levels[levels != 0],
                        extend='both', cmap='seismic', transform=ccrs.PlateCarree())
        plt.colorbar(cf, ax=ax, orientation="horizontal", pad=.05, extendrect=True, ticks=levels)
        plt.title(str(height)+"hPa, Omega (ps s-1)")
    else:
        var = var*1000.
        cf = ax.contourf(to_np(lons), to_np(lats), var, levels,
                        extend='both', cmap='seismic', transform=ccrs.PlateCarree())
        plt.colorbar(cf, ax=ax, orientation="horizontal", pad=.05, extendrect=True, ticks=levels)
        plt.title(str(height)+"hPa, Vertical wind speed (m h-1)") # (m s-1)")

    if message == None:
        message = var_name+"_"+str(height)+"hPa"
    else:
        message = message+"_"+var_name+"_"+str(height)+"hPa"

    fig.savefig('./plots/5Nov/wind_subsidence/3Nov/spatial_map_wrf_wind_subsidence_'+message , bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    #######################################################
    # Decks to run:
    #    heat_advection
    #######################################################
    case_name  = "hw2009_3Nov" # "hw2013_3Nov"
    var_name   = "wa" # "omg" #
    var_unit   = "km h-1" # None #
    height     = 850

    if case_name == "hw2009_3Nov":
        period     = "20090122-20090213"
        time_s = datetime(2009,1,28,0,0,0,0)
        time_e = datetime(2009,2,8,23,59,0,0)
    elif  case_name == "hw2013_3Nov":
        period     = "20121229-20130122"
        time_s = datetime(2013,1,4,0,0,0,0)
        time_e = datetime(2013,1,18,23,59,0,0)
    elif  case_name == "hw2019_3Nov":
        period     = "20190108-20190130"
        time_s = datetime(2019,1,14,0,0,0)
        time_e = datetime(2019,1,26,23,59,0,0)

    cpl_atmo_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_'+period+'_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_'+period+'_fd'  # atmo output of wrf-cable run

    file_paths        = [cpl_atmo_file_gw] #cpl_atmo_file_fd, cpl_atmo_file_fd,cpl_atmo_file_fd, cpl_atmo_file_gw

    if len(file_paths) > 1:
        message = "Couple_GW-FD_"+str(time_s)+"-"+str(time_e)
    else:
        message = "Couple_GW_"+str(time_s)+"-"+str(time_e)

    wind_subsidence(file_paths, var_name, var_unit, height, time_s, time_e, message=message) #  loc_lat=loc_lat, loc_lon=loc_lat,
