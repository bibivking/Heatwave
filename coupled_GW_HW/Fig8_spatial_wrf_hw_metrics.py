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

'''
copy from script_submission/spatial_wrf_hw_metrics.py but only use AWAP Tmax
'''

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

def hw_thresholds(AWAP_tmax_file, AWAP_tmax90_file, percent):

    method = "calc"

    if method == "calc":
        year_tot       = 2019 - 1970 + 1
        Time_tmp,tmax  = read_var(AWAP_tmax_file, "tmax", lat_name="latitude", lon_name="longitude")

        tmax_day    = np.delete(tmax,np.s_[365*2+31+28::365*3+366],0)
        tmax_day_4D = tmax_day.reshape(year_tot, 365, np.shape(tmax)[1], np.shape(tmax)[2])

        tmax_percent= np.percentile(tmax_day_4D,percent,axis=0)
        tmax        = None
        # ========== write out percentile ==========
        # read lat and lon
        AWAP_tmax       = nc.Dataset(AWAP_tmax_file, 'r')
        lat             = AWAP_tmax.variables['latitude'][:]
        lon             = AWAP_tmax.variables['longitude'][:]
        nday            = 365
        nlon            = len(lon)
        nlat            = len(lat)

        # create file and write global attributes
        f               = nc.Dataset(AWAP_tmax90_file, 'w', format='NETCDF4')
        f.description   = str(percent)+' th percentile of 1970-2019 daily max temperature, created by MU Mengyuan'
        f.source        = AWAP_tmax_file
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

        tmax            = f.createVariable('tmax', 'f4', ('time', 'latitude','longitude',))
        tmax.units     = "C"
        tmax.long_name = str(percent) + "th of air temperature"
        tmax.CF_name   = str(percent) + "th of air temperature"
        tmax.proj4     = "+proj=longlat +ellps=GRS80"
        tmax[:]        = tmax_percent

        # print("np.any(np.isnan(tmax))")
        # print(np.any(np.isnan(tmax)))
        # print(tmax)

        f.close()

    elif method == "plot":

        # check mask
        file = nc.Dataset(AWAP_tmax90_file, mode='r')
        lon  = file.variables['longitude'][:]
        lat  = file.variables['latitude'][:]
        tmax = file.variables['tmax'][:]

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

            var_contours = plt.contourf(lon, lat, tmax[i,:,:], levels=np.arange(20,45,1), transform=ccrs.PlateCarree(), cmap=cmap)
            plt.colorbar(var_contours, ax=ax1, orientation="horizontal", pad=.05)
            plt.savefig('./plots/tmax_'+str(percent)+'th_day-'+str(i+1)+'.png',dpi=300)

def regrid_AWAP_to_WRF(AWAP_file, WRF_file):

    awap_file = Dataset(AWAP_file, mode='r')
    tmax      = awap_file.variables["tmax"][:]
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
        value = np.reshape(tmax[day,:,:],-1)
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

def make_hw_days_file(AWAP_tmax_file, AWAP_tmax90_file, AWAP_masked_file, WRF_file):

    print("======== In mask_hw_days =========")
    # =============== Read AWAP tmax  ===============
    AWAP_max      = Dataset(AWAP_tmax_file, mode='r')
    time          = nc.num2date(AWAP_max.variables['time'][:],AWAP_max.variables['time'].units,
                    only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    tmax          = AWAP_max.variables['tmax'][:]
    lats          = AWAP_max.variables["latitude"][:]
    lons          = AWAP_max.variables["longitude"][:]
    nlat          = len(lats)
    nlon          = len(lons)
    ntime         = len(time)

    # =============== Read hw threshold ===============
    AWAP_tmax90   = Dataset(AWAP_tmax90_file, mode='r')
    # Time          = nc.num2date(AWAP_tmax90.variables['time'][:],AWAP_tmax90.variables['time'].units,
    #                 only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    tmax_90th     = AWAP_tmax90.variables["tmax"][:]
    # tmax_90th     = regrid_AWAP_to_WRF_var(Tmax_90th, AWAP_tmax90_file, WRF_file)

    # =============== Make new array ===============
    tmax_mask      = np.zeros([ntime - 12, nlat,nlon]) # minus 12 Feb 29

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
            # print("Time[doy] ",Time[doy])
            tmax_mask[cnt,:,:] = \
                np.where(tmax[i,:,:] >= tmax_90th[doy,:,:], tmax[i,:,:], np.nan)
            time_new.append(i)
            cnt   = cnt+1
        else:
            print(the_time.year, "is a leap year")
            if doy < 59:
                print("before Feb 29")
                print("doy=",doy)
                # print("Time[doy] ",Time[doy])
                # before Feb 29
                tmax_mask[cnt,:,:] = \
                    np.where(tmax[i,:,:] >= tmax_90th[doy,:,:], tmax[i,:,:], np.nan)
                time_new.append(i)
                cnt  = cnt+1
            elif doy > 59:
                print("after Feb 29")
                # after Feb 29
                doy = doy - 1
                print("doy=",doy)
                # print("Time[doy] ",Time[doy])
                tmax_mask[cnt,:,:] = \
                    np.where(tmax[i,:,:] >= tmax_90th[doy,:,:], tmax[i,:,:], np.nan)
                time_new.append(i)
                cnt  = cnt+1

    # =============== create the file ===============
    ntime = ntime - 12 # change to new ntime

    # regrid to AWAP
    tmax_mask_regrid     = regrid_AWAP_to_WRF_var(tmax_mask, AWAP_tmax90_file, WRF_file)
    wrf_file  = Dataset(WRF_file, mode='r')
    lat_out   = wrf_file.variables["XLAT"][0,:,:]
    lon_out   = wrf_file.variables["XLONG"][0,:,:]
    nlat      = len(lat_out[:,0])
    nlon      = len(lon_out[0,:])

    # create file and write global attributes
    f               = nc.Dataset(AWAP_masked_file, 'w', format='NETCDF4')
    f.description   = 'The days tmax temperature above ' \
                    + '90 th percentile of 1970-2019 daily max temperature. ' \
                    + 'The map has been resampled to the wrf domain. The dataset skipped each Feb 29. ' \
                    + 'Created by MU Mengyuan.'

    f.source        = AWAP_tmax90_file + ', ' + AWAP_tmax_file
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
    time.units      = AWAP_max.variables['time'].units
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
    latitude[:,:]   = lat_out

    longitude       = f.createVariable('longitude', 'f4', ('latitude','longitude'))
    longitude.units = "degrees_east"
    longitude.long_name = "longitude"
    longitude.standard_name = "projection_x_coordinate"
    longitude.axis  = "X"
    longitude[:,:]  = lon_out

    tmax           = f.createVariable('tmax', 'f4', ('time', 'latitude','longitude',))
    tmax.units     = "C"
    tmax.long_name = "air temperature above 90th daily max temperature"
    tmax.CF_name   = "air temperature above 90th daily max temperature"
    tmax.proj4     = "Lambert Conformal"
    tmax[:]        = tmax_mask_regrid

    f.close()

#
# def plot_spatial_T_abv_thrshld(AWAP_tmax90_file, file_path, month_s, tot_day=30, loc_lat=None, loc_lon=None,
#                                lat_names=None, lon_names=None, is_wrf=False, AWAP_mean_file=None):
#
#     print("======== In plot_spatial_T_abv_thrshld =========")
#     if is_wrf:
#         # =============== Read wrf data ================
#         # Open the NetCDF file
#         ncfile   = Dataset(file_path)
#         ntime    = len(ncfile.variables['Times'][:,0])
#
#         # read WRF time
#         time_tmp = []
#         for i in np.arange(ntime):
#             time_temp = datetime.strptime(str(ncfile.variables['Times'][i,:], 'utf-8'),'%Y-%m-%d_%H:%M:%S')
#             time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1))
#         time  = np.array(time_tmp)
#
#         # to get lat and lon
#         p     = getvar(ncfile, "pressure")
#
#         # Get the lat/lon coordinates
#         lats, lons = latlon_coords(p)
#
#         # Get the cartopy mapping object
#         cart_proj = get_cartopy(p)
#
#         # Extract the pressure, geopotential height, and wind variables
#         Var   = read_wrf_surf_var(file_path, 'T2', loc_lat, loc_lon)
#         Var   = Var - 273.15
#         # =============== Read hw threshold ===============
#         tmax   = regrid_AWAP_to_WRF(AWAP_tmax90_file, file_path)
#     else:
#         # =============== Read hw threshold ===============
#         if AWAP_mean_file == None:
#             awap_file  = Dataset(AWAP_tmax90_file, mode='r')
#             lat        = awap_file.variables["latitude"][:]
#             lon        = awap_file.variables["longitude"][:]
#             lons, lats = np.meshgrid(lon, lat)
#             tmax        = awap_file.variables["tmax"][:]
#         else:
#             WRF_file   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg/wrfout_20090122-20090213_gw"
#             wrf_file   = Dataset(WRF_file, mode='r')
#             lats       = wrf_file.variables["XLAT"][0,:,:]
#             lons       = wrf_file.variables["XLONG"][0,:,:]
#             tmax        = regrid_AWAP_to_WRF(AWAP_tmax90_file, WRF_file)
#
#             AWAP_90th  = Dataset(AWAP_tmax90_file, mode='r')
#             time_90th  = nc.num2date(AWAP_90th.variables['time'][:],AWAP_90th.variables['time'].units,
#                             only_use_cftime_datetimes=False, only_use_python_datetimes=True)
#
#         # =============== Read AWAP tmax ================
#         # Open the NetCDF file
#         if AWAP_mean_file == None:
#             # Calc AWAP tmax
#             AWAP_tmax_file = "/g/data/w97/mm3972/data/AWAP/AWAP_AUS_temp/AWAP_daily_tmax_1970_2019.nc"
#             AWAP_tmin_file = "/g/data/w97/mm3972/data/AWAP/AWAP_AUS_temp/AWAP_daily_tmin_1970_2019.nc"
#
#             time, tmax    = read_var(AWAP_tmax_file, "tmax", loc_lat, loc_lon, "latitude", "longitude")
#             time, tmin    = read_var(AWAP_tmin_file, "tmin", loc_lat, loc_lon, "latitude", "longitude")
#             Var           = (tmax+tmin)/2.
#         else:
#             # Read AWAP tmax
#             AWAP_mean     = Dataset(AWAP_mean_file, mode='r')
#             time          = nc.num2date(AWAP_mean.variables['time'][:],AWAP_mean.variables['time'].units,
#                             only_use_cftime_datetimes=False, only_use_python_datetimes=True)
#             time          = time - datetime(2000,1,1,0,0,0)
#             Var           = AWAP_mean.variables['tmax'][:]
#
#     # =============== calc T > threshold ===============
#     for i in np.arange(tot_day):
#         print("In calc T > threshold")
#         time_s  = month_s + timedelta(days=int(i))
#         time_e  = month_s + timedelta(days=int(i+1)) - timedelta(seconds=1)
#         doy     = time_s  - datetime(time_s.year,1,1,0,0,0,0)
#         print("time_s",time_s)
#         print("time_e",time_e)
#         print("int(doy.days)",int(doy.days))
#         print("time_90th[doy.days]",time_90th[doy.days])
#
#         var     = spital_var(time, Var, time_s, time_e)
#         var_tmp = np.where(var >= tmax[int(doy.days),:,:], var, -999.)
#
#         # ============= Plot setting ===============
#         # Create the figure
#         fig = plt.figure(figsize=(12,9))
#         fig.subplots_adjust(hspace=0.3)
#         fig.subplots_adjust(wspace=0.2)
#
#         plt.rcParams['text.usetex']     = False
#         plt.rcParams['font.family']     = "sans-serif"
#         plt.rcParams['font.serif']      = "Helvetica"
#         plt.rcParams['axes.linewidth']  = 1.5
#         plt.rcParams['axes.labelsize']  = 14
#         plt.rcParams['font.size']       = 14
#         plt.rcParams['legend.fontsize'] = 12
#         plt.rcParams['xtick.labelsize'] = 12
#         plt.rcParams['ytick.labelsize'] = 14
#
#         almost_black = '#262626'
#         # change the tick colors also to the almost black
#         plt.rcParams['ytick.color']     = almost_black
#         plt.rcParams['xtick.color']     = almost_black
#
#         # change the text colors also to the almost black
#         plt.rcParams['text.color']      = almost_black
#
#         # Change the default axis colors from black to a slightly lighter black,
#         # and a little thinner (0.5 instead of 1)
#         plt.rcParams['axes.edgecolor']  = almost_black
#         plt.rcParams['axes.labelcolor'] = almost_black
#
#         # set the box type of sequence number
#         props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
#         # choose colormap
#         if is_wrf:
#             # Set the GeoAxes to the projection used by WRF
#             ax = plt.axes(projection=cart_proj)
#         else:
#             ax = plt.axes(projection=ccrs.PlateCarree())
#         # Download and add the states and coastlines
#         states = NaturalEarthFeature(category="cultural", scale="50m",
#                                              facecolor="none",
#                                              name="admin_1_states_provinces_shp")
#         ax.add_feature(states, linewidth=.5, edgecolor="black")
#         ax.coastlines('50m', linewidth=0.8)
#
#         # start plotting
#         if loc_lat == None:
#             ax.set_extent([135,155,-40,-25])
#         else:
#             ax.set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])
#
#         # Add the var contours
#         levels = [16.,18.,20.,22.,24.,26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.]
#
#         var_contours = plt.contourf(to_np(lons), to_np(lats), to_np(var_tmp),
#                                     transform=ccrs.PlateCarree(), cmap=get_cmap("seismic"),
#                                     levels = levels, extend='both') #,"jet" #“rainbow”#"coolwarm"
#         plt.colorbar(var_contours, ax=ax, orientation="horizontal", pad=.05)  #[levels!=0]"bwr" [levels!=0] levels = levels[levels!=0],
#
#         if is_wrf:
#             # Set the map bounds
#             ax.set_xlim(cartopy_xlim(p))
#             ax.set_ylim(cartopy_ylim(p))
#
#         plt.title(str(time_s))
#         if is_wrf:
#             fig.savefig('./plots/figures/HW_cover/spatial_map_tmax90th_HW_cover_WRF_'+str(time_s)+'.png', bbox_inches='tight', pad_inches=0.1)
#         else:
#             fig.savefig('./plots/figures/HW_cover/spatial_map_tmax90th_HW_cover_AWAP_'+str(time_s)+'_read_tmax.png', bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    # ###################################
    # Deck make_hw_days_file
    # ###################################

    AWAP_tmax_file     = "/g/data/w97/mm3972/data/AWAP/AWAP_AUS_temp/AWAP_daily_tmax_1970_2019.nc"
    AWAP_tmax90_file   = "/g/data/w97/mm3972/scripts/Heatwave/coupled_GW_HW/nc_file/AWAP_tmax_1970_2019_90-th.nc"
    AWAP_masked_file   = "/g/data/w97/mm3972/scripts/Heatwave/coupled_GW_HW/nc_file/AWAP_tmax_above_tmax90th_1970-2019_WRF_domain.nc"
    WRF_file           = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg/wrfout_20090122-20090213_gw"
    # hw_thresholds(AWAP_tmax_file, AWAP_tmax90_file, 90)
    make_hw_days_file(AWAP_tmax_file, AWAP_tmax90_file, AWAP_masked_file, WRF_file)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
