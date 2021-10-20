#!/usr/bin/python

'''
Functions:
1. Analyze weather situation during heatwaves
2. Process ERAI, AWAP, offline CABLE and LIS-CABLE data
'''

from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from wrf import (getvar, interplevel, get_cartopy, cartopy_xlim,
                 cartopy_ylim, to_np, latlon_coords)
from common_utils import *

def plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=None, loc_lon=None, lat_names=None, lon_names=None, message=None):

    print("======== In plot_spital_map =========")

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    time1, Var1  = read_var(file_paths[0], var_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
    time1, lats1 = read_var(file_paths[0], lat_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
    time1, lons1 = read_var(file_paths[0], lon_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])

    if var_names[0] in ['tas','Tair','Tair_f_inst']:
        var1         = spital_var(time1,Var1,time_s,time_e)-273.15
    elif var_names[0] in ['tp']:
        scale        = get_scale(var_names[0])
        var1         = spital_ERAI_tp(time1,Var1,time_s,time_e)*scale
    else:
        scale        = get_scale(var_names[0])
        var1         = spital_var(time1,Var1,time_s,time_e)*scale

    if len(file_paths) > 1:
        time2, Var2  = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
        time2, lats2 = read_var(file_paths[1], lat_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
        time2, lons2 = read_var(file_paths[1], lon_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
        scale        = get_scale(var_names[1])
        var2         = spital_var(time2,Var2,time_s,time_e)*scale

    if len(file_paths) > 2:
        time3, Var3  = read_var(file_paths[2], var_names[2], loc_lat, loc_lon, lat_names[2], lon_names[2])
        time3, lats3 = read_var(file_paths[2], lat_names[2], loc_lat, loc_lon, lat_names[2], lon_names[2])
        time3, lons3 = read_var(file_paths[2], lon_names[2], loc_lat, loc_lon, lat_names[2], lon_names[2])
        scale        = get_scale(var_names[2])
        var3         = spital_var(time3,Var3,time_s,time_e)*scale

    if len(file_paths) > 3:
        time4, Var4  = read_var(file_paths[3], var_names[3], loc_lat, loc_lon, lat_names[3], lon_names[3])
        time4, lats4 = read_var(file_paths[3], lat_names[3], loc_lat, loc_lon, lat_names[3], lon_names[3])
        time4, lons4 = read_var(file_paths[3], lon_names[3], loc_lat, loc_lon, lat_names[3], lon_names[3])
        scale        = get_scale(var_names[3])
        var4         = spital_var(time4,Var4,time_s,time_e)*scale

    fig = plt.figure(figsize=(6,5))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # start plotting
    if loc_lat == None:
        # ax.set_extent([140,154,-40,-28])
        ax.set_extent([135,155,-40,-25])
    else:
        ax.set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])

    ax.coastlines(resolution="50m",linewidth=1)

    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
    gl.xlabels_top   = False
    gl.ylabels_right = False
    gl.xlines        = True

    if loc_lat == None:
        # gl.xlocator = mticker.FixedLocator([140,145,150])
        # gl.ylocator = mticker.FixedLocator([-40,-35,-30])
        gl.xlocator     = mticker.FixedLocator([135,140,145,150,155])
        gl.ylocator     = mticker.FixedLocator([-40,-35,-30,-25])
    else:
        gl.xlocator = mticker.FixedLocator(loc_lon)
        gl.ylocator = mticker.FixedLocator(loc_lat)

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':'black'}
    gl.ylabel_style = {'size':10, 'color':'black'}


    # plot Var1
    if var_names[0] in ['tas','Tair','Tair_f_inst']:
        clevs = np.linspace( 15.,45., num=31)
        cmap  = plt.cm.hot_r
    else:
        # clevs = np.linspace( 0.,120., num=13)
        clevs = np.linspace( 0.,20., num=5)
        cmap  = plt.cm.GnBu # BrBG

    plt.contourf(lons1, lats1, var1, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #,#bwr)#coolwarm)#cm.BrBG) # clevs,

    plt.title(var_names[0], size=16)
    cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    # cb.set_label(units,size=14,rotation=270,labelpad=15)
    cb.ax.tick_params(labelsize=10)

    # plot Var2
    if len(file_paths) > 1 and var_names[1] == 'ps':
        clevs = np.linspace( 900.,1100., num=11)
        cs = plt.contour(lons2, lats2, var2, clevs, transform=ccrs.PlateCarree(),linewidths=0.8,colors="darkgray") #, ,cmap=plt.cm.hot_r)#bwr)#coolwarm)#cm.BrBG) # clevs,
        cl = plt.clabel(cs, inline=True, fmt="%4d",fontsize=6) #manual=True)

    # plot Var3, Var4
    if len(file_paths) > 3 and var_names[2] == 'uas' and var_names[3] == 'vas':
        qv = plt.quiver(lons1[::3,::3], lats1[::3,::3], var3[::3,::3], var4[::3,::3], scale=200, color='k')

    if message == None:
        message = var_names[0]
    else:
        message = message + "_" + var_names[0]
    plt.savefig('./plots/19Oct/weather/spatial_map_weather_analysis_'+message+'.png',dpi=300)

def plot_spatial_map_hgt(file_paths, var_name, var_unit, height, time_s, time_e, loc_lat=None, loc_lon=None, message=None):

    # Open the NetCDF file
    encoding = 'utf-8' # Times in WRF output is btype, convert to string

    ncfile1  = Dataset(file_paths[0])
    ntime    = len(ncfile1.variables['Times'][:,0])
    time_tmp = []
    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(ncfile1.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(time_temp - datetime(2000,1,1))
    print(type(time_tmp))
    time = np.array(time_tmp)
    print(time)
    # Extract the pressure, geopotential height, and wind variables
    p1    = getvar(ncfile1, "pressure")
    Var1  = read_wrf_var(file_paths[0], var_name, var_unit, height, loc_lat, loc_lon)
    Z1    = read_wrf_var(file_paths[0], "z", "dm", height, loc_lat, loc_lon)
    Ua1   = read_wrf_var(file_paths[0], "ua", "m s-1", height, loc_lat, loc_lon)
    Va1   = read_wrf_var(file_paths[0], "va", "m s-1", height, loc_lat, loc_lon)

    if var_name in ['temp']:
        var1  = spital_var(time,Var1,time_s,time_e)-273.15
    else:
        scale = get_scale(var_name)
        var1  = spital_var(time,Var1,time_s,time_e)*scale

    z1   = spital_var(time,Z1,time_s,time_e)
    ua1  = spital_var(time,Ua1,time_s,time_e)
    va1  = spital_var(time,Va1,time_s,time_e)

    if len(file_paths) > 1:
        ncfile2 = Dataset(file_paths[1])

        p2    = getvar(ncfile1, "pressure")
        var2  = read_wrf_var(file_paths[1], var_name, var_unit, height, loc_lat, loc_lon)
        z2    = read_wrf_var(file_paths[1], "z", "dm", height, loc_lat, loc_lon)
        ua2   = read_wrf_var(file_paths[1], "ua", "m s-1", height, loc_lat, loc_lon)
        va2   = read_wrf_var(file_paths[1], "va", "m s-1", height, loc_lat, loc_lon)

        if var_name in ['temp']:
            var2  = spital_var(time,Var2,time_s,time_e)-273.15
        else:
            scale = get_scale(var_name)
            var2  = spital_var(time,Var2,time_s,time_e)*scale

        z2   = spital_var(time,Z2,time_s,time_e)
        ua2  = spital_var(time,Ua2,time_s,time_e)
        va2  = spital_var(time,Va2,time_s,time_e)

        # Calculate difference
        var_diff = var_hgt2 - var_hgt1
        u_diff   = u_hgt2 - u_hgt1
        v_diff   = v_hgt2 - v_hgt1
        z_diff   = z_hgt2 - z_hgt1

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

    # Convert the lat/lon coordinates to x/y coordinates in the projection space
    # x, y = bm(to_np(lons), to_np(lats))

    # Add the 500 hPa geopotential height contours
    # if height == 500:
    #     levels = np.arange(-10., 10, 2.)
    # elif height == 850:
    # levels = np.arange(100., 10., 2.)
    contours = plt.contour(to_np(lons), to_np(lats), to_np(z1),  colors="black", transform=ccrs.PlateCarree())#levels = levels,
    plt.clabel(contours, inline=1, fontsize=6, fmt="%i")

    # Add the var contours
    # levels = np.arange(val_min, val_max, 0.1)
    var_contours = plt.contourf(to_np(lons), to_np(lats), to_np(var1),
                         transform=ccrs.PlateCarree(), cmap=get_cmap("bwr"),extend='both') #levels = levels,"jet" #“rainbow”#"coolwarm"
    plt.colorbar(var_contours, ax=ax, orientation="horizontal", pad=.05)

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(p1))
    ax.set_ylim(cartopy_ylim(p1))

    # Add the 500 hPa wind barbs, only plotting every 125th data point.
    ax.quiver(to_np(lons[::5,::5]), to_np(lats[::5,::5]), to_np(ua1[::5, ::5]),
             to_np(va1[::5, ::5]), scale=20., transform=ccrs.PlateCarree()) # width=0.0002,

    plt.title(str(height)+"hPa Geopotential Height (dm), Temperature (degC) and Barbs (m s-1)")

    if message == None:
        message = var_name+'_'+str(height)+"hPa"
    else:
        message = message+"_"+var_name+'_'+str(height)+"hPa"

    fig.savefig('./plots/19Oct/weather/spatial_map_weather_analysis_wrf_'+message , bbox_inches='tight', pad_inches=0.1)


if __name__ == "__main__":


    # ======================= Option =======================
    region = "SE Aus"

    # ====================== Pre-load =======================
    ERAI_path    = '/g/data/ub4/erai/netcdf/3hr/atmos/oper_fc_sfc/v01'
    ERAI_T_file  = ERAI_path + '/tas/tas_3hrs_ERAI_historical_fc-sfc_20190101_20190131.nc' # air temperature
    ERAI_P_file  = ERAI_path + '/ps/ps_3hrs_ERAI_historical_fc-sfc_20190101_20190131.nc'   # surface pressure
    ERAI_U_file  = ERAI_path + '/uas/uas_3hrs_ERAI_historical_fc-sfc_20190101_20190131.nc' # 10 m wind speed
    ERAI_V_file  = ERAI_path + '/vas/vas_3hrs_ERAI_historical_fc-sfc_20190101_20190131.nc' # 10 m wind speed
    ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20190101_20190131.nc' # Total rainfall

    AWAP_path    = '/g/data/w35/Shared_data/AWAP_3h_v1'
    AWAP_T_file  = AWAP_path + '/Tair/AWAP.Tair.3hr.2019.nc'  # air temperature
    AWAP_R_file  = AWAP_path + '/Rainf/AWAP.Rainf.3hr.2019.nc'  # Daily rainfall

    cpl_land_path  = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/gw_jan2019/LIS_output'
    cpl_land_file  = cpl_land_path + '/LIS.CABLE.201901-201901.d01.nc'  # land output of wrf-cable run

    cpl_atmo_file = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/gw_jan2019/WRF_output'
    cpl_atmo_file  = cpl_atmo_file + '/wrfout_d01_2019-01-01_00:00:00'  # atmo output of wrf-cable run

    # datetime(year, month, day, hour, minute, second, microsecond)
    hw_s    = [ datetime(2009,1,28,0,0,0,0),
                datetime(2011,1,30,0,0,0,0),
                datetime(2013,1,1,0,0,0,0),
                datetime(2019,1,12,0,0,0,0)
                ]
                # datetime(2014,1,13,0,0,0,0),
                # datetime(2017,1,31,0,0,0,0),
    hw_e    = [ datetime(2009,2,8,0,0,0,0),
                datetime(2011,2,6,0,0,0,0),
                datetime(2013,1,9,0,0,0,0),
                datetime(2019,1,26,0,0,0,0)
                ]
                # datetime(2014,1,18,0,0,0,0),
                # datetime(2017,2,12,0,0,0,0),
    if region == "Aus":
        loc_lat    = [-44,-10]
        loc_lon    = [112,154]
    elif region == "SE Aus":
        loc_lat    = [-40,-25]
        loc_lon    = [135,155]
    # =================== Opreation ====================

    # '''
    # Plot ERAI air temp + 10 m wind speed + pressure
    # '''
    # file_paths  = [ERAI_T_file, ERAI_P_file, ERAI_U_file, ERAI_V_file]

    # var_names   = ['tas','ps','uas','vas']
    # lat_names   = ["lat","lat","lat","lat"]
    # lon_names   = ["lon","lon","lon","lon"]

    # for i in np.arange(0,30):
    #     time_s = datetime(2019,1,1,0,0,0,0) + timedelta(days=int(i))
    #     time_e = datetime(2019,1,2,0,0,0,0) + timedelta(days=int(i))
    #     message = "ERAI_2019-01-"+str(i+1)
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)

    #################################
    #    Plot ERAI rainfall         #
    #################################

    # # 2009
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20090101_20090131.nc' # Total rainfall
    # time_s = datetime(2009,1,1,0,0,0,0)
    # time_e = datetime(2009,2,1,0,0,0,0)
    # message = "ERAI_2009-1_Rainfall"

    # # 2011
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20110101_20110131.nc' # Total rainfall
    # time_s = datetime(2011,1,1,0,0,0,0)
    # time_e = datetime(2011,2,1,0,0,0,0)
    # message = "ERAI_2011-1_Rainfall"

    # # 2013
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20121201_20121231.nc' # Total rainfall
    # time_s = datetime(2012,12,1,0,0,0,0)
    # time_e = datetime(2013,1,1,0,0,0,0)
    # message = "ERAI_2012-12_Rainfall"

    # # 2019
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20181201_20181231.nc' # Total rainfall
    # time_s = datetime(2018,12,1,0,0,0,0)
    # time_e = datetime(2019,1,1,0,0,0,0)
    # message = "ERAI_2018-12_Rainfall"

    hw_e    = [ datetime(2009,2,8,0,0,0,0),
                datetime(2011,2,6,0,0,0,0),
                datetime(2013,1,9,0,0,0,0),
                datetime(2019,1,26,0,0,0,0)
                ]


    # # 2009
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20090201_20090228.nc' # Total rainfall
    # time_s = datetime(2009,2,1,0,0,0,0)
    # time_e = hw_e[0]
    # message = "ERAI_2009hw_Rainfall"

    # # 2011
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20110201_20110228.nc' # Total rainfall
    # time_s = datetime(2011,2,1,0,0,0,0)
    # time_e = hw_e[1]
    # message = "ERAI_2011hw_Rainfall"

    # # 2013
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20130101_20130131.nc' # Total rainfall
    # time_s = datetime(2013,1,1,0,0,0,0)
    # time_e = hw_e[2]
    # message = "ERAI_2013hw_Rainfall"

    # 2019
    ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20190101_20190131.nc' # Total rainfall
    time_s = datetime(2019,1,1,0,0,0,0)
    time_e = datetime(2019,1,12,0,0,0,0)#hw_e[3]
    message = "ERAI_20190101-12_Rainfall"

    file_paths  = [ERAI_R_file]

    var_names   = ['tp']
    lat_names   = ["lat"]
    lon_names   = ["lon"]

    # for i in np.arange(0,30):
    #     time_s = datetime(2019,1,1,0,0,0,0) + timedelta(days=int(i))
    #     time_e = datetime(2019,1,2,0,0,0,0) + timedelta(days=int(i))
    #     message = "ERAI_2019-01-"+str(i+1)
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)

    plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
                    lon_names=lon_names,message=message)

    # '''
    # Plot AWAP Tair
    # '''

    # file_paths  = [AWAP_T_file]

    # var_names   = ['Tair']
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]

    # for i in np.arange(0,30):
    #     time_s = datetime(2019,1,1,0,0,0,0) + timedelta(days=int(i))
    #     time_e = datetime(2019,1,2,0,0,0,0) + timedelta(days=int(i))
    #     message = "AWAP_2019-01-"+str(i+1)
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)

    # '''
    # Plot AWAP rainfall
    # '''

    # file_paths  = [AWAP_R_file]

    # var_names   = ['Rainf']
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]

    # for i in np.arange(0,30):
    #     time_s = datetime(2019,1,1,0,0,0,0) + timedelta(days=int(i))
    #     time_e = datetime(2019,1,2,0,0,0,0) + timedelta(days=int(i))
    #     message = "AWAP_2019-01-"+str(i+1)
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)


    # '''
    # Plot WRF-CABLE land var
    # '''

    # file_paths  = [cpl_land_file]

    # var_names   = ['Tair_f_inst']#'Rainf_tavg','Tair_f_inst','Qair_f_inst','Psurf_f_inst']
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]

    # for i in np.arange(0,30):
    #     time_s = datetime(2019,1,1,0,0,0,0) + timedelta(days=int(i))
    #     time_e = datetime(2019,1,2,0,0,0,0) + timedelta(days=int(i))
    #     message = "Couple_2019-01-"+str(i+1)
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)

    # '''
    # Plot WRF-CABLE atmo var
    # '''
    #
    # file_paths  = [cpl_atmo_file]
    #
    # var_name    = "temp"
    # var_unit    = "degC"
    # height      = 500
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]
    #
    # for i in np.arange(0,30):
    #     time_s = datetime(2019,1,1,0,0,0,0) + timedelta(days=int(i))
    #     time_e = datetime(2019,1,2,0,0,0,0) + timedelta(days=int(i))
    #     message = "Couple_2019-01-"+str(i+1)
    #     plot_spatial_map_hgt(file_paths, var_name, var_unit, height, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lat, message=message)
