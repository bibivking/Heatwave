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
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                        cartopy_ylim, latlon_coords, ALL_TIMES)
from common_utils import *

def get_time_cood(file_path, Time_s ,Time_e):

    # ****************** process time ******************
    ncfile = Dataset(file_path)

    Time   = nc.num2date(ncfile.variables['time'][:],ncfile.variables['time'].units,
                only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time   = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)

    time_s = Time_s - datetime(2000,1,1,0,0,0)
    time_e = Time_e - datetime(2000,1,1,0,0,0)

    time_cood = (time>=time_s) & (time<time_e)

    return time_cood

def hw_thresholds(AWAP_tmax_file, AWAP_tmin_file, AWAP_tmn90_file, percent):

    method = "plot"

    if method == "calc":
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
        f               = nc.Dataset(AWAP_tmn90_file, 'w', format='NETCDF4')
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

    elif method == "plot":

        # check mask
        file = nc.Dataset(AWAP_tmn90_file, mode='r')
        lon  = file.variables['longitude'][:]
        lat  = file.variables['latitude'][:]
        tmn  = file.variables['tmean'][:]

        cmap = plt.cm.coolwarm # #,"jet" #“rainbow”#"coolwarm"

        for i in np.arange(365):
            print(i, "-day")

            fig1 = plt.figure(figsize=(12,9))
            # Get the lat/lon coordinates

            ax1 = plt.axes(projection=ccrs.PlateCarree())
            ax1.set_extent([110,155,-45,-10])
            ax1.coastlines(resolution="50m",linewidth=1)

            # Add gridlines
            gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
            gl.xlabels_top  = False
            gl.ylabels_right= False
            gl.xlines       = True
            gl.xlocator     = mticker.FixedLocator([110,115,120,125,130,135,140,145,150,155])
            gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15,-10])
            gl.xformatter   = LONGITUDE_FORMATTER
            gl.yformatter   = LATITUDE_FORMATTER
            gl.xlabel_style = {'size':10, 'color':'black','rotation': 90}
            gl.ylabel_style = {'size':10, 'color':'black'}

            var_contours = plt.contourf(lon, lat, tmn[i,:,:], levels=np.arange(20,45,1), transform=ccrs.PlateCarree(), cmap=cmap)
            plt.colorbar(var_contours, ax=ax1, orientation="horizontal", pad=.05)
            plt.savefig('./plots/5Nov/hw_metric/3Nov/Tmean_'+str(percent)+'th_day-'+str(i+1)+'.png',dpi=300)

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

def regrid_AWAP_to_WRF_var(AWAP_var, AWAP_file, WRF_file):
    
    awap_file = Dataset(AWAP_file, mode='r')
    lat_in    = awap_file.variables["latitude"][:]
    lon_in    = awap_file.variables["longitude"][:]

    wrf_file  = Dataset(WRF_file, mode='r')
    lat_out   = wrf_file.variables["XLAT"][0,:,:]
    lon_out   = wrf_file.variables["XLONG"][0,:,:]
    nlat      = len(lat_out[:,0])
    nlon      = len(lon_out[0,:])
    print(nlat)
    print(nlon)

    lon_in_2D, lat_in_2D = np.meshgrid(lon_in,lat_in)
    lon_in_1D            = np.reshape(lon_in_2D,-1)
    lat_in_1D            = np.reshape(lat_in_2D,-1)
    
    ntime                = np.shape(AWAP_var[:,0,0])[0]
    print("ntime = ",ntime)
    
    Value                = np.zeros([ntime,nlat,nlon])

    for day in np.arange(ntime):
        value = np.reshape(AWAP_var[day,:,:],-1)
        Value[day,:,:] =  griddata((lon_in_1D, lat_in_1D), value, (lon_out, lat_out), method="nearest")

    print("np.shape(Value)")
    print(np.shape(Value))

    return Value

def make_hw_days_file(AWAP_tmax_file, AWAP_tmin_file, AWAP_tmn90_file, AWAP_mean_file, 
                      AWAP_masked_file, WRF_file, make_Tmn=True):
    
    print("======== In mask_hw_days =========")
    
    if make_Tmn:
        
        # =============== Read WRF lat and lon dimensions ===============
        wrf_file  = Dataset(WRF_file, mode='r')
        lat_out   = wrf_file.variables["XLAT"][0,:,:]
        lon_out   = wrf_file.variables["XLONG"][0,:,:]
        nlat      = len(lat_out[:,0])
        nlon      = len(lon_out[0,:])
        
        # =============== Read AWAP time dimension tmean =================
        # Open the NetCDF file
        time, tmax    = read_var(AWAP_tmax_file, "tmax", lat_name="latitude", lon_name="longitude")
        time, tmin    = read_var(AWAP_tmin_file, "tmin", lat_name="latitude", lon_name="longitude")
        Var           = (tmax+tmin)/2.
        Tmn           = regrid_AWAP_to_WRF_var(Var, AWAP_tmn90_file, WRF_file)
        
        AWAP_file     = Dataset(AWAP_tmax_file, mode='r')
        time          = nc.num2date(AWAP_file.variables['time'][:],AWAP_file.variables['time'].units,
                        only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        ntime         = len(time)
        
        # =============== create the file ===============
        
        # create file and write global attributes
        f             = nc.Dataset(AWAP_mean_file, 'w', format='NETCDF4')
        f.description = 'AWAP 1970-2019 daily mean temperature regridded to the WRF domain. Created by MU Mengyuan.'
                            
        f.source        = AWAP_tmax_file + " and " + AWAP_tmin_file
        f.history       = "Created by: %s" % (os.path.basename(__file__))
        f.creation_date = "%s" % (datetime.now())

        # set dimensions
        f.createDimension('time', ntime)
        f.createDimension('latitude', nlat)
        f.createDimension('longitude', nlon)
        f.Conventions   = "CF-1.0"

        # create variables
        time            = f.createVariable('time', 'f4', ('time',))
        time.units      = AWAP_file.variables['time'].units
        time.long_name  = "time"
        time.standard_name = "time"
        time.calendar   = "standard"
        time.axis       = "T"
        time[:]         = np.arange(ntime)

        latitude        = f.createVariable('latitude', 'f4', ('latitude','longitude'))
        latitude.units  = "degrees_north"
        latitude.long_name = "latitude"
        latitude.standard_name = "projection_y_coordinate"
        latitude.axis   = "Y"
        latitude[:,:]   = lat_out

        longitude       = f.createVariable('longitude', 'f4', ('latitude','longitude'))
        longitude.units = "degrees_east"
        longitude.long_name = "longitude"
        longitude.standard_name = "projection_x_coordinate"
        longitude.axis  = "X"
        longitude[:,:]  = lon_out

        tmean           = f.createVariable('tmean', 'f4', ('time', 'latitude','longitude',))
        tmean.units     = "C"
        tmean.long_name = "daily mean temperature"
        tmean.CF_name   = "daily mean temperature"
        tmean.proj4     = "Lambert Conformal"
        tmean[:]        = Tmn

        f.close()
    
        print("check whether time read by read_var is corret: ",time)
        
    else:
        # =============== Read AWAP Tmean  ===============
        AWAP_mean     = Dataset(AWAP_mean_file, mode='r')
        time          = nc.num2date(AWAP_mean.variables['time'][:],AWAP_mean.variables['time'].units,
                        only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        Tmn           = AWAP_mean.variables['tmean'][:]
        lats          = AWAP_mean.variables["latitude"][:,:]
        lons          = AWAP_mean.variables["longitude"][:,:]
        nlat          = len(lats[:,0])
        nlon          = len(lons[0,:])
        ntime         = len(time)
        
        # =============== Read hw threshold ===============
        awap_file     = Dataset(AWAP_tmn90_file, mode='r')
        Time          = nc.num2date(awap_file.variables['time'][:],awap_file.variables['time'].units,
                        only_use_cftime_datetimes=False, only_use_python_datetimes=True) 
        tmean         = awap_file.variables["tmean"][:]
        Tmn_90th      = regrid_AWAP_to_WRF_var(tmean, AWAP_tmn90_file, WRF_file)

        # =============== Make new array ===============
        Tmn_mask      = np.zeros([ntime - 12, nlat,nlon]) # minus 12 Feb 29
        
        # =============== calc T > threshold ===============
        cnt           = 0 
        time_new      = []
        for i in np.arange(ntime):
            print("time[i] ",time[i])
            the_time  = time[i]
            delta_day = the_time - datetime(the_time.year,1,1,0,0,0,0)
            doy       = int(delta_day.days)
            print("year = ",the_time.year)
            print("day of year = ",doy)
            if (the_time.year % 4) != 0:
                print(the_time.year, "isn't a leap year")
                print("Time[doy] ",Time[doy])
                Tmn_mask[cnt,:,:] = \
                    np.where(Tmn[i,:,:] >= Tmn_90th[doy,:,:], Tmn[i,:,:], np.nan)
                time_new.append(i)
                cnt   = cnt+1
            else: 
                print(the_time.year, "is a leap year")
                if doy < 59:
                    print("before Feb 29")
                    print("doy=",doy)
                    print("Time[doy] ",Time[doy])
                    # before Feb 29
                    Tmn_mask[cnt,:,:] = \
                        np.where(Tmn[i,:,:] >= Tmn_90th[doy,:,:], Tmn[i,:,:], np.nan)
                    time_new.append(i)
                    cnt  = cnt+1
                elif doy > 59:
                    print("after Feb 29")
                    # after Feb 29
                    doy = doy - 1
                    print("doy=",doy)
                    print("Time[doy] ",Time[doy])
                    Tmn_mask[cnt,:,:] = \
                        np.where(Tmn[i,:,:] >= Tmn_90th[doy,:,:], Tmn[i,:,:], np.nan)
                    time_new.append(i)
                    cnt  = cnt+1        
                            
        # # =============== calc T > threshold ===============
        # for i in np.arange(tot_day):

        #     time_s  = month_s + timedelta(days=int(i))
        #     time_e  = month_s + timedelta(days=int(i+1)) - timedelta(seconds=1)
        #     doy     = time_s  - datetime(time_s.year,1,1,0,0,0,0)

        #     var     = spital_var(time, Var, time_s, time_e)
        #     var_tmp = np.where(var >= Tmn[int(doy.days),:,:], var, -999.)
            
        # =============== create the file ===============
        ntime = ntime - 12 # change to new ntime
        
        # create file and write global attributes
        f               = nc.Dataset(AWAP_masked_file, 'w', format='NETCDF4')
        f.description   = 'The days mean temperature above ' \
                        + '90 th percentile of 1970-2019 daily mean temperature. ' \
                        + 'The map has been resampled to the wrf domain. The dataset skipped each Feb 29. ' \
                        + 'Created by MU Mengyuan.'
                            
        f.source        = AWAP_tmn90_file + ', ' + AWAP_mean_file
        f.history       = "Created by: %s" % (os.path.basename(__file__))
        f.creation_date = "%s" % (datetime.now())

        # set dimensions
        f.createDimension('time', ntime)
        f.createDimension('latitude', nlat)
        f.createDimension('longitude', nlon)
        f.Conventions  = "CF-1.0"

        print(time_new)
        # create variables
        time            = f.createVariable('time', 'f4', ('time',))
        time.units      = AWAP_mean.variables['time'].units
        time.long_name  = "time"
        time.standard_name = "time"
        time.calendar   = "standard"
        time.axis       = "T"
        time[:]         = time_new

        latitude        = f.createVariable('latitude', 'f4', ('latitude','longitude'))
        latitude.units  = "degrees_north"
        latitude.long_name = "latitude"
        latitude.standard_name = "projection_y_coordinate"
        latitude.axis   = "Y"
        latitude[:,:]   = lats

        longitude       = f.createVariable('longitude', 'f4', ('latitude','longitude'))
        longitude.units = "degrees_east"
        longitude.long_name = "longitude"
        longitude.standard_name = "projection_x_coordinate"
        longitude.axis  = "X"
        longitude[:,:]  = lons

        tmean           = f.createVariable('tmean', 'f4', ('time', 'latitude','longitude',))
        tmean.units     = "C"
        tmean.long_name = "air temperature above 90th daily mean temperature"
        tmean.CF_name   = "air temperature above 90th daily mean temperature"
        tmean.proj4     = "Lambert Conformal"
        tmean[:]        = Tmn_mask

        f.close()

def count_days_abv_thrshld(var, thrshld, time, time_s, time_e):

    time_cood = time_mask(time, time_s, time_e)

    time_tmp  = time + datetime(2000,1,1)
    doy       = time_tmp - datetime(time_tmp[-1].year,1,1,0,0,0,0)

    hours     = np.zeros(np.shape(var[0,:,:]))
    print(np.shape(var))

    for i in np.arange(len(time_cood)):
        if time_cood[i]:
            hours = hours + np.where( var[int(i),:,:]- 273.15 >= thrshld[int(doy[i].days),:,:], 1, 0)

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
    p    = getvar(ncfile1, "pressure", timeidx=ALL_TIMES)

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
    lats, lons = latlon_coords(p)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(p)

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
                   levels = levels, transform=ccrs.PlateCarree(), cmap=get_cmap("seismic"),extend='both') #,"jet" #“rainbow”#"coolwarm"
    plt.colorbar(var_contours, ax=ax, orientation="horizontal", pad=.05)  #"bwr"

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(p))
    ax.set_ylim(cartopy_ylim(p))

    plt.title(metric)

    if message == None:
        message = var_name+"_"+metric
    else:
        message = message+"_"+var_name+"_"+metric

    fig.savefig('./plots/5Nov/hw_metric/3Nov/spatial_wrf_hw_metrics_'+message , bbox_inches='tight', pad_inches=0.1)

def plot_spatial_wrf_hw_durition(AWAP_tmax_file, AWAP_tmin_file, AWAP_tmn90_file, file_paths, var_name, time_s,
                                 time_e, percent=90, metric='Tmean',loc_lat=None, loc_lon=None, message=None):

    ### Calculate AWAP heatwave threshold ###
    # hw_thresholds(AWAP_tmax_file, AWAP_tmin_file, AWAP_tmn90_file, percent)

    # ### regrid AWAP to WRF domain ###
    Tmn_regrid = regrid_AWAP_to_WRF(AWAP_tmn90_file, file_paths[0])

    # Open the NetCDF file
    ncfile1  = Dataset(file_paths[0])
    ntime    = len(ncfile1.variables['Times'][:,0])
    time_tmp = []

    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(ncfile1.variables['Times'][i,:], 'utf-8'),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1))

    time  = np.array(time_tmp)

    # to get lat and lon
    p    = getvar(ncfile1, "pressure", timeidx=ALL_TIMES)

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
    lats, lons = latlon_coords(p)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(p)

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
        # levels = [-3.9,-2.9,-1.9,-0.9,0.9,1.9,2.9,3.9]
        # np.arange(-10., 11.,2.)
    else:
        levels = np.arange(0., 50.,1.)
    print(levels)

    var_contours = plt.contourf(to_np(lons), to_np(lats), to_np(days),
                                transform=ccrs.PlateCarree(), cmap=get_cmap("seismic"),
                                levels = levels, extend='both') #,"jet" #“rainbow”#"coolwarm"
    plt.colorbar(var_contours, ax=ax, orientation="horizontal", pad=.05)  #[levels!=0]"bwr" [levels!=0] levels = levels[levels!=0],
    #
    # Set the map bounds
    ax.set_xlim(cartopy_xlim(p))
    ax.set_ylim(cartopy_ylim(p))

    plt.title(metric)

    if message == None:
        message = var_name+"_"+metric
    else:
        message = message+"_"+var_name+"_"+metric

    fig.savefig('./plots/5Nov/hw_metric/3Nov/spatial_wrf_hw_duriation_'+message , bbox_inches='tight', pad_inches=0.1)

def plot_spatial_T_abv_thrshld(AWAP_tmn90_file, file_path, month_s, tot_day=30, loc_lat=None, loc_lon=None, 
                               lat_names=None, lon_names=None, is_wrf=False, AWAP_mean_file=None):

    print("======== In plot_spatial_T_abv_thrshld =========")
    if is_wrf:
        # =============== Read wrf data ================
        # Open the NetCDF file
        ncfile   = Dataset(file_path)
        ntime    = len(ncfile.variables['Times'][:,0])

        # read WRF time
        time_tmp = []
        for i in np.arange(ntime):
            time_temp = datetime.strptime(str(ncfile.variables['Times'][i,:], 'utf-8'),'%Y-%m-%d_%H:%M:%S')
            time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1))
        time  = np.array(time_tmp)

        # to get lat and lon
        p     = getvar(ncfile, "pressure")

        # Get the lat/lon coordinates
        lats, lons = latlon_coords(p)

        # Get the cartopy mapping object
        cart_proj = get_cartopy(p)

        # Extract the pressure, geopotential height, and wind variables
        Var   = read_wrf_surf_var(file_path, 'T2', loc_lat, loc_lon)
        Var   = Var - 273.15
        # =============== Read hw threshold ===============
        Tmn   = regrid_AWAP_to_WRF(AWAP_tmn90_file, file_path)
    else:
        # =============== Read hw threshold ===============
        if AWAP_mean_file == None:
            awap_file  = Dataset(AWAP_tmn90_file, mode='r')
            lat        = awap_file.variables["latitude"][:]
            lon        = awap_file.variables["longitude"][:]
            lons, lats = np.meshgrid(lon, lat)
            Tmn        = awap_file.variables["tmean"][:]
        else:
            WRF_file   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg/wrfout_20090122-20090213_gw"
            wrf_file   = Dataset(WRF_file, mode='r')
            lats       = wrf_file.variables["XLAT"][0,:,:]
            lons       = wrf_file.variables["XLONG"][0,:,:]
            Tmn        = regrid_AWAP_to_WRF(AWAP_tmn90_file, WRF_file)
            
            AWAP_90th  = Dataset(AWAP_tmn90_file, mode='r')
            time_90th  = nc.num2date(AWAP_90th.variables['time'][:],AWAP_90th.variables['time'].units,
                            only_use_cftime_datetimes=False, only_use_python_datetimes=True)
   
        # =============== Read AWAP Tmean ================
        # Open the NetCDF file
        if AWAP_mean_file == None:
            # Calc AWAP Tmean
            AWAP_tmax_file = "/g/data/w35/mm3972/data/AWAP/AWAP_AUS_temp/AWAP_daily_tmax_1970_2019.nc"
            AWAP_tmin_file = "/g/data/w35/mm3972/data/AWAP/AWAP_AUS_temp/AWAP_daily_tmin_1970_2019.nc"

            time, tmax    = read_var(AWAP_tmax_file, "tmax", loc_lat, loc_lon, "latitude", "longitude")
            time, tmin    = read_var(AWAP_tmin_file, "tmin", loc_lat, loc_lon, "latitude", "longitude")
            Var           = (tmax+tmin)/2.
        else:
            # Read AWAP Tmean
            AWAP_mean     = Dataset(AWAP_mean_file, mode='r')
            time          = nc.num2date(AWAP_mean.variables['time'][:],AWAP_mean.variables['time'].units,
                            only_use_cftime_datetimes=False, only_use_python_datetimes=True)
            time          = time - datetime(2000,1,1,0,0,0)
            Var           = AWAP_mean.variables['tmean'][:]
            
    # =============== calc T > threshold ===============
    for i in np.arange(tot_day):
        print("In calc T > threshold")
        time_s  = month_s + timedelta(days=int(i))
        time_e  = month_s + timedelta(days=int(i+1)) - timedelta(seconds=1)
        doy     = time_s  - datetime(time_s.year,1,1,0,0,0,0)
        print("time_s",time_s)
        print("time_e",time_e)
        print("int(doy.days)",int(doy.days))
        print("time_90th[doy.days]",time_90th[doy.days])

        var     = spital_var(time, Var, time_s, time_e)
        var_tmp = np.where(var >= Tmn[int(doy.days),:,:], var, -999.)

        # ============= Plot setting ===============
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
        if is_wrf:
            # Set the GeoAxes to the projection used by WRF
            ax = plt.axes(projection=cart_proj)
        else:
            ax = plt.axes(projection=ccrs.PlateCarree())
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

        # Add the var contours
        levels = [16.,18.,20.,22.,24.,26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.]

        var_contours = plt.contourf(to_np(lons), to_np(lats), to_np(var_tmp),
                                    transform=ccrs.PlateCarree(), cmap=get_cmap("seismic"),
                                    levels = levels, extend='both') #,"jet" #“rainbow”#"coolwarm"
        plt.colorbar(var_contours, ax=ax, orientation="horizontal", pad=.05)  #[levels!=0]"bwr" [levels!=0] levels = levels[levels!=0],

        if is_wrf:
            # Set the map bounds
            ax.set_xlim(cartopy_xlim(p))
            ax.set_ylim(cartopy_ylim(p))

        plt.title(str(time_s))
        if is_wrf:
            fig.savefig('./plots/figures/HW_cover/spatial_map_Tmn90th_HW_cover_WRF_'+str(time_s)+'.png', bbox_inches='tight', pad_inches=0.1)
        else:
            fig.savefig('./plots/figures/HW_cover/spatial_map_Tmn90th_HW_cover_AWAP_'+str(time_s)+'_read_Tmn.png', bbox_inches='tight', pad_inches=0.1)

def plot_spatial_EHF(file_path, time_s, time_e, loc_lat=None, loc_lon=None, lat_names=None, lon_names=None):

    print("======== In plot_spital_map =========")

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    time_cood   = get_time_cood(file_path, time_s ,time_e)
    print(time_cood)
    Time, Event = read_var(file_path, "event", loc_lat, loc_lon, "lat", "lon")
    print("np.shape(Event)",np.shape(Event))
    print("Time",Time)

    time  = Time[time_cood] + datetime(2000,1,1,0,0,0)

    event = Event[time_cood]
    print("time",time)
    print("np.shape(event)",np.shape(event))

    ncfile      = Dataset(file_path)
    lat         = ncfile.variables['lat'][:]
    lon         = ncfile.variables['lon'][:]
    print("lat",lat)
    print("lon",lon)

    for i in np.arange(len(time)):

        fig = plt.figure(figsize=(6,5))
        ax  = plt.axes(projection=ccrs.PlateCarree())

        # start plotting
        if loc_lat == None:
            ax.set_extent([130,155,-45,-20])
        else:
            ax.set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])

        ax.coastlines(resolution="50m",linewidth=1)

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top   = False
        gl.ylabels_right = False
        gl.xlines        = True

        if loc_lat == None:
            gl.xlocator     = mticker.FixedLocator([135,140,145,150,155])
            gl.ylocator     = mticker.FixedLocator([-40,-35,-30,-25])
        else:
            gl.xlocator = mticker.FixedLocator(loc_lon)
            gl.ylocator = mticker.FixedLocator(loc_lat)

        gl.xformatter   = LONGITUDE_FORMATTER
        gl.yformatter   = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black'}
        gl.ylabel_style = {'size':10, 'color':'black'}

        # cmap  = plt.cm.hot_r
        levels= [0.1,0.9]
        plt.contourf(lon, lat, event[i], levels=levels, transform=ccrs.PlateCarree(), extend='both') # cmap=cmap,
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

        cb.ax.tick_params(labelsize=10)

        plt.savefig('./plots/figures/spatial_map_EHF_HW_cover_'+str(time[i])+'.png',dpi=300)
        fig = None
        ax  = None

if __name__ == "__main__":
    
    # ###################################
    # Deck make_hw_days_file
    # ###################################
    
    AWAP_tmax_file     = "/g/data/w35/mm3972/data/AWAP/AWAP_AUS_temp/AWAP_daily_tmax_1970_2019.nc"
    AWAP_tmin_file     = "/g/data/w35/mm3972/data/AWAP/AWAP_AUS_temp/AWAP_daily_tmin_1970_2019.nc"
    AWAP_tmn90_file    = "/g/data/w35/mm3972/scripts/Heatwave/coupled_GW_HW/nc_file/AWAP_tmean_1970_2019_90-th.nc"
    AWAP_mean_file     = "/g/data/w35/mm3972/scripts/Heatwave/coupled_GW_HW/nc_file/AWAP_tmean_1970-2019_WRF_domain.nc"
    AWAP_masked_file   = "/g/data/w35/mm3972/scripts/Heatwave/coupled_GW_HW/nc_file/AWAP_tmean_above_Tmn90th_1970-2019_WRF_domain.nc"
    WRF_file           = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg/wrfout_20090122-20090213_gw"
    make_hw_days_file(AWAP_tmax_file, AWAP_tmin_file, AWAP_tmn90_file, AWAP_mean_file, AWAP_masked_file, WRF_file, make_Tmn=False)

    # # #######################################
    # # Deck for other than make_hw_days_file
    # # #######################################
    # var_name     = 'T2'
    # case_names   = ["hw2019_3Nov"]# "hw2009_3Nov",  "hw2013_3Nov", 

    # for case_name in case_names:
    #     if case_name == "hw2009_3Nov":
    #         # 2009 HW : 28–31.01.2009 and 6–8.02.2009
    #         period = "20090122-20090213"
    #         Time_s = datetime(2009,1,28,0,0,0,0)
    #         Time_e = datetime(2009,2,8,23,59,0,0)
    #         month_s= datetime(2009,1,1,0,0,0,0)
    #         # Time_s = datetime(2009,1,22,0,0,0,0)
    #         # Time_e = datetime(2009,2,13,23,59,0,0)
    #     elif  case_name == "hw2013_3Nov":
    #         # 2013 HW : 4–8.01.2013, ​11–13.01.2013 and 17–18.01.2013
    #         period = "20121229-20130122"
    #         Time_s = datetime(2013,1,4,0,0,0,0)
    #         Time_e = datetime(2013,1,18,23,59,0,0)
    #         month_s= datetime(2013,1,1,0,0,0,0)
    #         # Time_s = datetime(2012,12,29,0,0,0,0)
    #         # Time_e = datetime(2013,1,22,23,59,0,0)
    #     elif  case_name == "hw2019_3Nov":
    #         # 2019 HW : 14–18.01.2019 and 22–26.01.2019
    #         period = "20190108-20190130"
    #         Time_s = datetime(2019,1,14,0,0,0)
    #         Time_e = datetime(2019,1,26,23,59,0,0)
    #         month_s= datetime(2019,1,1,0,0,0,0)
    #         # Time_s = datetime(2019,1,8,0,0,0)
    #         # Time_e = datetime(2019,1,30,23,59,0,0)

    #     #######################################################
    #     # Decks to run:
    #     #    plot_spatial_EHF
    #     #######################################################
    #     AWAP_tmn90_file = "./nc_file/AWAP_tmean_1970_2019_90-th.nc"
    #     AWAP_mean_file  = "/g/data/w35/mm3972/scripts/Heatwave/coupled_GW_HW/nc_file/AWAP_tmean_1970-2019_WRF_domain.nc"
    #     cpl_atmo_file   = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'
    #     file_path       = cpl_atmo_file + '/wrfout_'+period+'_gw'  # atmo output of wrf-cable run
    #     tot_day         = 61
    #     plot_spatial_T_abv_thrshld(AWAP_tmn90_file, file_path, month_s, tot_day=tot_day, AWAP_mean_file=AWAP_mean_file)
        
    #     # ######################################################
    #     # Decks to run:
    #     #    plot_spatial_EHF
    #     # ######################################################
    #     # EHF_file = '/g/data/w35/mm3972/scripts/ehfheatwaves/nc_file/AUS_based_on_tmean_95th/EHF_heatwaves_daily_2009-2019.nc'
        
    #     # plot_spatial_EHF(EHF_file, Time_s, Time_e)


    #     #######################################################
    #     # Decks to run:
    #     #    plot_spatial_wrf_hw_durition
    #     #######################################################

    #     # AWAP_tmax_file = "/g/data/w35/mm3972/data/AWAP/AWAP_AUS_temp/AWAP_daily_tmax_1970_2019.nc"
    #     # AWAP_tmin_file = "/g/data/w35/mm3972/data/AWAP/AWAP_AUS_temp/AWAP_daily_tmin_1970_2019.nc"

    #     # cpl_atmo_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'
    #     # cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_'+period+'_gw'  # atmo output of wrf-cable run
    #     # cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_'+period+'_fd'  # atmo output of wrf-cable run

    #     # file_paths        = [ cpl_atmo_file_fd, cpl_atmo_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw

    #     # percent       = 90

    #     # AWAP_tmn90_file = "./nc_file/AWAP_tmean_1970_2019_"+str(percent)+"-th.nc"

    #     # tot_day   = (Time_e-Time_s).days + 1

    #     # # for i in np.arange(tot_day):
    #     #     # time_s = Time_s + timedelta(days=int(i))
    #     #     # time_e = Time_s + timedelta(days=int(i+1)) - timedelta(seconds=1)

    #     # if len(file_paths) > 1:
    #     #     message = 'Couple_GW-FD_'+str(Time_s)+'-'+str(Time_e)
    #     # else:
    #     #     message = 'Couple_GW_'+str(Time_s)+'-'+str(Time_e)

    #     # plot_spatial_wrf_hw_durition(AWAP_tmax_file, AWAP_tmin_file, AWAP_tmn90_file, file_paths, var_name, Time_s, Time_e, percent=percent, message=message)

    #     #######################################################
    #     # Decks to run:
    #     #    plot_spatial_wrf_hw_metrics
    #     #######################################################

    #     # if len(file_paths) > 1:
    #     #     message = 'Couple_GW-FD_'+str(time_s)+'-'+str(time_e)
    #     # else:
    #     #     message = 'Couple_GW_'+str(time_s)+'-'+str(time_e)

    #     # metric   = 'Tmean'
    #     # plot_spatial_wrf_hw_metrics(file_paths, var_name, time_s, time_e, metric=metric, message=message)

    #     # metric   = 'Tmax'
    #     # plot_spatial_wrf_hw_metrics(file_paths, var_name, time_s, time_e, metric=metric, message=message)
