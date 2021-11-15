#!/usr/bin/python

import os
import sys
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.animation as animation
from scipy.interpolate import griddata
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                        cartopy_ylim, latlon_coords, ALL_TIMES)
from common_utils import *

def hw_thresholds(AWAP_tmax_file, AWAP_tmin_file, AWAP_file_out, percent):

    year_tot       = 2019 - 1970 + 1
    Time_tmp,tmax  = read_var(AWAP_tmax_file, "tmax", lat_name="latitude", lon_name="longitude")
    Time_tmp,tmin  = read_var(AWAP_tmin_file, "tmin", lat_name="latitude", lon_name="longitude")

    tmn        = (tmax + tmin)/2.

    tmn_day    = np.delete(tmn,np.s_[365*2+31+28::365*3+366],0)

    tmn_day_4D = tmn_day.reshape(year_tot, 365, np.shape(tmn)[1], np.shape(tmn)[2])

    tmn_percent= np.percentile(tmn_day_4D,percent,axis=0)

    # ========== write out percentile ==========
    # read lat and lon
    AWAP_tmax       = nc.Dataset(AWAP_tmax_file, 'r')
    lat             = AWAP_tmax.variables['latitude'][:]
    lon             = AWAP_tmax.variables['longitude'][:]
    nday            = 365
    nlon            = len(lon)
    nlat            = len(lat)

    # create file and write global attributes
    f               = nc.Dataset(AWAP_file_out, 'w', format='NETCDF4')
    f.description   = str(percent)+' th percentile of 1970-2019 daily mean temperature, created by MU Mengyuan'
    f.source        = AWAP_tmax_file + " and " + AWAP_tmin_file
    f.history       = "Created by: %s" % (os.path.basename(__file__))
    f.creation_date = "%s" % (datetime.now())

    # set dimensions
    f.createDimension('time', nday)
    f.createDimension('latitude', nlat)
    f.createDimension('longitude', nlon)
    f.Conventions  = "CF-1.0"

    # create variables
    time            = f.createVariable('time', 'f4', ('time',))
    time.units      = "days since 2019-01-01 00:00:00"
    time.long_name  = "time"
    time.standard_name = "time"
    time.calendar   = "standard"
    time.axis       = "T"
    time[:]         = np.arange(365)
    # print(time)

    latitude        = f.createVariable('latitude', 'f4', ('latitude'))
    latitude.units  = "degrees_north"
    latitude.long_name = "latitude"
    latitude.standard_name = "projection_y_coordinate"
    latitude.axis   = "Y"
    latitude[:]     = lat
    # print(latitude)

    # print("np.any(np.isnan(latitude))")
    # print(np.any(np.isnan(latitude)))

    longitude       = f.createVariable('longitude', 'f4', ('longitude'))
    longitude.units = "degrees_east"
    longitude.long_name = "longitude"
    longitude.standard_name = "projection_x_coordinate"
    longitude.axis  = "X"
    longitude[:]    = lon
    # print(longitude)
    #
    # print("np.any(np.isnan(longitude))")
    # print(np.any(np.isnan(longitude)))

    tmean           = f.createVariable('tmean', 'f4', ('time', 'latitude','longitude',))
    tmean.units     = "C"
    tmean.long_name = str(percent) + "th of air temperature"
    tmean.CF_name   = str(percent) + "th of air temperature"
    tmean.proj4     = "+proj=longlat +ellps=GRS80"
    tmean[:]        = tmn_percent

    # print("np.any(np.isnan(tmean))")
    # print(np.any(np.isnan(tmean)))
    # print(tmean)

    f.close()

def regrid_AWAP_to_WRF(AWAP_file, WRF_file):

    awap_file = Dataset(AWAP_file, mode='r')
    tmean     = awap_file.variables["tmean"][:]
    lat_in    = awap_file.variables["latitude"][:]
    lon_in    = awap_file.variables["longitude"][:]

    wrf_file  = Dataset(WRF_file, mode='r')
    lat_out   = wrf_file.variables["XLAT"][0,:,:]
    lon_out   = wrf_file.variables["XLONG"][0,:,:]
    nlat      = len(lat_out[:,0])
    nlon      = len(lon_out[0,:])

    lon_in_2D, lat_in_2D = np.meshgrid(lon_in,lat_in)
    lon_in_1D            = np.reshape(lon_in_2D,-1)
    lat_in_1D            = np.reshape(lat_in_2D,-1)

    Value                = np.zeros((365,nlat,nlon))
    print(np.shape(Value))

    for day in np.arange(365):
        value = np.reshape(tmean[day,:,:],-1)
        Value[day,:,:] =  griddata((lon_in_1D, lat_in_1D), value, (lon_out, lat_out), method="nearest")

    print("np.shape(Value)")
    print(np.shape(Value))

    return Value

def count_days_abv_thrshld(var, thrshld, time, time_s, time_e):

    time_cood = time_mask(time, time_s, time_e)

    time_tmp  = time + datetime(2000,1,1)
    doy       = time_tmp - datetime(time_tmp[-1].year,1,1,0,0,0,0)

    hours     = np.zeros(np.shape(var[0,:,:]))
    print(np.shape(var))

    for i in np.arange(len(time_cood)):
        if time_cood[i]:
            hours = hours + np.where( var[int(i),:,:]- 273.15 >= thrshld[int(doy[i].days),:,:], 1, 0)

    print(hours[4,6])

    return hours

def plot_spatial_wrf_hw_metrics(file_paths, var_name, time_s, time_e, metric='Tmean',loc_lat=None, loc_lon=None, message=None):

    # Open the NetCDF file
    encoding = 'utf-8' # Times in WRF output is btype, convert to string

    ncfile1  = Dataset(file_paths[0])
    ntime    = len(ncfile1.variables['Times'][:,0])
    time_tmp = []

    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(ncfile1.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1))

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
            levels = [-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3,3.5,4]
            #np.arange(-3.0, 3.0, 0.2)
        elif metric == 'Tmax':
            levels = [-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3,3.5,4]
    else:
        levels = np.arange(np.nanmin(var), np.nanmax(var), 21)

    var_contours = plt.contourf(to_np(lons), to_np(lats), to_np(var),
                   levels = levels, transform=crs.PlateCarree(), cmap=get_cmap("seismic"),extend='both') #,"jet" #“rainbow”#"coolwarm"
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

def plot_spatial_wrf_hw_durition(AWAP_tmax_file, AWAP_tmin_file, AWAP_file_out, file_paths, var_name, time_s, time_e, percent=90, metric='Tmean',loc_lat=None, loc_lon=None, message=None):

    ### Calculate AWAP heatwave threshold ###
    # hw_thresholds(AWAP_tmax_file, AWAP_tmin_file, AWAP_file_out, percent)

    ### regrid AWAP to WRF domain ###
    Tmn_regrid = regrid_AWAP_to_WRF(AWAP_file_out, file_paths[0])

    # Open the NetCDF file
    ncfile1  = Dataset(file_paths[0])
    ntime    = len(ncfile1.variables['Times'][:,0])
    time_tmp = []

    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(ncfile1.variables['Times'][i,:], 'utf-8'),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1))

    time  = np.array(time_tmp)

    # to get lat and lon
    p1    = getvar(ncfile1, "pressure", timeidx=ALL_TIMES)

    # Extract the pressure, geopotential height, and wind variables
    Var1  = read_wrf_surf_var(file_paths[0], var_name, loc_lat, loc_lon)
    Days1 = count_days_abv_thrshld(Var1, Tmn_regrid, time, time_s, time_e)

    if len(file_paths) > 1:
        Var2  = read_wrf_surf_var(file_paths[1], var_name, loc_lat, loc_lon)
        Days2 = count_days_abv_thrshld(Var2, Tmn_regrid, time, time_s, time_e)
        days  = Days2 - Days1
    else:
        days  = Days1

    print(days)

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
        levels = [-20.,-18.,-16.,-14.,-12.,-10.,-8.,-6.,-4.,-2.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20.]
        # np.arange(-10., 11.,2.)
    else:
        levels = np.arange(0., 50.,1.)
    print(levels)

    var_contours = plt.contourf(to_np(lons), to_np(lats), to_np(days),
                                transform=crs.PlateCarree(), cmap=get_cmap("seismic"),
                                levels = levels, extend='both') #,"jet" #“rainbow”#"coolwarm"
    plt.colorbar(var_contours, ax=ax, orientation="horizontal", pad=.05)  #[levels!=0]"bwr" [levels!=0] levels = levels[levels!=0],
    #
    # Set the map bounds
    ax.set_xlim(cartopy_xlim(p1))
    ax.set_ylim(cartopy_ylim(p1))

    plt.title(metric)

    if message == None:
        message = var_name+"_"+metric
    else:
        message = message+"_"+var_name+"_"+metric

    fig.savefig('./plots/5Nov/hw_metric/3Nov/spatial_wrf_hw_duriation_'+message , bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    var_name       = 'T2'
    case_name      = "hw2019_3Nov" # "hw2013_3Nov"# "hw2019_3Nov"#

    if case_name == "hw2009_3Nov":
        # 2009 HW : 28–31.01.2009 and 6–8.02.2009
        period = "20090122-20090213"
        time_s = datetime(2009,1,28,0,0,0,0)
        time_e = datetime(2009,2,8,23,59,0,0)
    elif  case_name == "hw2013_3Nov":
        # 2013 HW : 4–8.01.2013, ​11–13.01.2013 and 17–18.01.2013
        period = "20121229-20130122"
        time_s = datetime(2013,1,4,0,0,0,0)
        time_e = datetime(2013,1,18,23,59,0,0)
    elif  case_name == "hw2019_3Nov":
        # 2019 HW : 14–18.01.2019 and 22–26.01.2019
        period = "20190108-20190130"
        time_s = datetime(2019,1,14,0,0,0)
        time_e = datetime(2019,1,26,23,59,0,0)

    AWAP_tmax_file = "/g/data/w35/mm3972/data/AWAP/AWAP_AUS_temp/AWAP_daily_tmax_1970_2019.nc"
    AWAP_tmin_file = "/g/data/w35/mm3972/data/AWAP/AWAP_AUS_temp/AWAP_daily_tmin_1970_2019.nc"

    cpl_atmo_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_'+period+'_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_'+period+'_fd'  # atmo output of wrf-cable run

    file_paths       = [cpl_atmo_file_fd, cpl_atmo_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw

    #######################################################
    # Decks to run:
    #    plot_spatial_wrf_hw_durition
    #######################################################
    percent       = 90

    AWAP_file_out = "./nc_file/AWAP_tmean_1970_2019_"+str(percent)+"-th.nc"

    if len(file_paths) > 1:
        message = 'Couple_GW-FD_'+str(time_s)+'-'+str(time_e)
    else:
        message = 'Couple_GW_'+str(time_s)+'-'+str(time_e)

    plot_spatial_wrf_hw_durition(AWAP_tmax_file, AWAP_tmin_file, AWAP_file_out, file_paths, var_name, time_s, time_e, percent=percent, message=message)

    #######################################################
    # Decks to run:
    #    plot_spatial_wrf_hw_metrics
    #######################################################

    # if len(file_paths) > 1:
    #     message = 'Couple_GW-FD_'+str(time_s)+'-'+str(time_e)
    # else:
    #     message = 'Couple_GW_'+str(time_s)+'-'+str(time_e)
    #
    # metric   = 'Tmean'
    # plot_spatial_wrf_hw_metrics(file_paths, var_name, time_s, time_e, metric=metric, message=message)
    #
    # metric   = 'Tmax'
    # plot_spatial_wrf_hw_metrics(file_paths, var_name, time_s, time_e, metric=metric, message=message)
