#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

def get_var_scale(var_name):

    '''
    Convert units
    '''

    var_s2d        = ["Rainf_tavg","Evap_tavg","ECanop_tavg","TVeg_tavg","ESoil_tavg"]
    var_umol_s2g_d = "GPP_tavg"

    s2d            = 3600*24. # s-1 to d-1
    umol_s2g_d     = 0.000001*12*s2d # umol s-1 to g d-1

    if var_name in var_s2d:
        scale = s2d
        units = "mm d-1"
    elif var_name in var_umol_s2g_d:
        scale = umol_s2g_d
        units = "gC m-2 d-1"
    else:
        scale = 1.
        units = None

    return (scale, units)

def plot_map_var(file_path, var_name):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var = Dataset(file_path, mode='r')

    time = var.variables['time']
    lons = var.variables['longitude'] # or ['lon']
    lats = var.variables['latitude']  # or ['lat']
    lon, lat = np.meshgrid(lons, lats)

    # 2-meter air temperatur K
    Var = var.variables[var_name]

    print(Var)

    # # Replace _FillValues with NaNs:
    # U10M_nans = U10M[:]
    # V10M_nans = V10M[:]
    # _FillValueU10M = U10M._FillValue
    # _FillValueV10M = V10M._FillValue
    # U10M_nans[U10M_nans == _FillValueU10M] = np.nan
    # V10M_nans[V10M_nans == _FillValueV10M] = np.nan

    # # Calculate wind speed:
    # ws = np.sqrt(U10M_nans**2+V10M_nans**2)
    # # Calculate wind direction in radians:
    # ws_direction = np.arctan2(V10M_nans,U10M_nans)

    # for i in np.arange(31):
        # NOTE: the MERRA-2 file contains hourly data for 24 hours (t=24). To get the daily mean wind speed, take the average of the hourly wind speeds:
        # T2M_daily_avg = np.nanmean(T2M[i*8:(i+1)*8,:,:], axis=0) -273.15 # i*8:(i+1)*8
        # print(T2M_daily_avg)
        # # NOTE: To calculate the average wind direction correctly it is important to use the 'vector average' as atan2(<v>,<u>) where <v> and <u> are the daily average component vectors, rather than as mean of the individual wind vector direction angle.  This avoids a situation where averaging 1 and 359 = 180 rather than the desired 0.
        # U10M_daily_avg = np.nanmean(U10M[i*8:(i+1)*8,:,:], axis=0)
        # V10M_daily_avg = np.nanmean(V10M[i*8:(i+1)*8,:,:], axis=0)

        # ws_daily_avg_direction = np.arctan2(V10M_daily_avg, U10M_daily_avg)

        #Plot Global MERRA-2 Wind Speed
        # Set the figure size, projection, and extent
        # fig = plt.figure(figsize=(8,4))
        # ax = plt.axes(projection=ccrs.Robinson())
        # ax.set_global()
        # ax.coastlines(resolution="110m",linewidth=1)
        # ax.gridlines(linestyle='--',color='black')
        # # Plot windspeed: set contour levels, then draw the filled contours and a colorbar
        # clevs = np.arange(0,19,1)
        # plt.contourf(lon, lat, ws_daily_avg, clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.jet)
        #
        # plt.title('MERRA-2 Daily Average 2-meter Wind Speed, 1 June 2010', size=14)
        # cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        # cb.set_label('m/s',size=12,rotation=0,labelpad=15)
        # cb.ax.tick_params(labelsize=10)
        # plt.savefig('Tair.png',dpi=1200)
        # plt.show()

    for i in np.arange(15675,15737):
        # The filled contours show the wind speed. The "quiver" function is used to overlay arrows to show the wind direction. The length of the arrows is determined by the wind speed.
        # Set the figure size, projection, and extent
        fig = plt.figure(figsize=(9,5))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([110,155,-45,-10])
        ax.coastlines(resolution="50m",linewidth=1)
        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator([110,115,120,125,130,135,140,145,150,155])
        gl.ylocator = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15,-10])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black'}
        gl.ylabel_style = {'size':10, 'color':'black'}
        # Plot windspeed
        clevs = np.arange(10,50,2)
        plt.contourf(lon, lat, Var[i,:,:], clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.jet) # T2M_daily_avg
        plt.title('AWAP Tmax ' + str(i-15675), size=16)
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.set_label('C',size=14,rotation=0,labelpad=15)
        cb.ax.tick_params(labelsize=10)
        # Overlay wind vectors
        # qv = plt.quiver(lon, lat, U10M[i,:,:], V10M[i,:,:], scale=420, color='k')
        plt.savefig('./plots/'+var_name+'-'+str(i-15675)+'.png',dpi=1200)

def plot_map_var_ts(file_path, wrf_path, var_names, ts):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var = Dataset(file_path, mode='r')
    wrf = Dataset(wrf_path,  mode='r')

    # use WRF output's lat & lon, since LIS output has default value
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    for var_name in var_names:

        print(var_name)
        scale,units = get_var_scale(var_name)
        Var = var.variables[var_name][ts]*scale

        # Make plots
        fig = plt.figure(figsize=(9,5))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([110,155,-45,-10])
        ax.coastlines(resolution="50m",linewidth=1)

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator([110,115,120,125,130,135,140,145,150,155])
        gl.ylocator = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15,-10])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black'}
        gl.ylabel_style = {'size':10, 'color':'black'}

        # clevs = np.arange(np.min(Var),np.max(Var)+(np.max(Var)-np.max(Var))/20.,(np.max(Var)-np.max(Var))/20.)
        if len(np.shape(Var)) == 2:
            plt.contourf(lon, lat, Var[:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.jet) # , T2M_daily_avg
            print("hi hi")
        if len(np.shape(Var)) == 3:
            plt.contourf(lon, lat, Var[0,:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.jet) # clevs, T2M_daily_avg

        plt.title(var_name+' ('+units+') ts='+str(ts), size=16)
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.set_label(units,size=14,rotation=0,labelpad=15)
        cb.ax.tick_params(labelsize=10)

        plt.savefig('./plots/'+var_name+'_ts-'+str(ts)+'.png',dpi=1200)
        Var = None

def plot_map_var_ts(file_path, wrf_path, var_names, ts):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var = Dataset(file_path, mode='r')
    wrf = Dataset(wrf_path,  mode='r')

    # use WRF output's lat & lon, since LIS output has default value
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    for var_name in var_names:

        print(var_name)
        scale, units = get_var_scale(var_name)
        Var = var.variables[var_name][ts]*scale

        # Make plots
        fig = plt.figure(figsize=(9,5))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([110,155,-45,-10])
        ax.coastlines(resolution="50m",linewidth=1)

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top  = False
        gl.ylabels_right= False
        gl.xlines       = True
        gl.xlocator     = mticker.FixedLocator([110,115,120,125,130,135,140,145,150,155])
        gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15,-10])
        gl.xformatter   = LONGITUDE_FORMATTER
        gl.yformatter   = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black'}
        gl.ylabel_style = {'size':10, 'color':'black'}

        # clevs = np.arange(np.min(Var),np.max(Var)+(np.max(Var)-np.max(Var))/20.,(np.max(Var)-np.max(Var))/20.)
        if len(np.shape(Var)) == 2:
            plt.contourf(lon, lat, Var[:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.jet) # , T2M_daily_avg
            print("hi hi")
        if len(np.shape(Var)) == 3:
            plt.contourf(lon, lat, Var[0,:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.jet) # clevs, T2M_daily_avg

        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        if units == None:
            plt.title(var_name+' ('+var.variables[var_name].units+') ts='+str(ts), size=16)
            cb.set_label(var.variables[var_name].units,size=14,rotation=0,labelpad=15)
        else:
            plt.title(var_name+' ('+units+') ts='+str(ts), size=16)
            cb.set_label(units,size=14,rotation=0,labelpad=15)
        cb.ax.tick_params(labelsize=10)

        plt.savefig('./plots/'+var_name+'_ts-'+str(ts)+'.png',dpi=1200)
        Var = None

def plot_map_var_diff_period_mean(file_paths, wrf_path, case_names, var_names, ts_s, ts_e):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var1 = Dataset(file_paths[0], mode='r')
    var2 = Dataset(file_paths[1], mode='r')

    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    for var_name in var_names:

        print(var_name)
        scale, units = get_var_scale(var_name)

        Var   = ( np.mean(var2.variables[var_name][ts_s:ts_e],axis=0)
                - np.mean(var1.variables[var_name][ts_s:ts_e],axis=0))*scale

        # Make plots
        fig = plt.figure(figsize=(9,5))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([110,155,-45,-10])
        ax.coastlines(resolution="50m",linewidth=1)

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator([110,115,120,125,130,135,140,145,150,155])
        gl.ylocator = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15,-10])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black'}
        gl.ylabel_style = {'size':10, 'color':'black'}

        # clevs = np.arange(np.min(Var),np.max(Var)+(np.max(Var)-np.max(Var))/20.,(np.max(Var)-np.max(Var))/20.)
        if len(np.shape(Var)) == 2:
            plt.contourf(lon, lat, Var[:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.seismic) # , T2M_daily_avg
            print("hi hi")
        if len(np.shape(Var)) == 3:
            plt.contourf(lon, lat, Var[0,:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.seismic) # clevs, T2M_daily_avg

        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        if units == None:
            plt.title(var_name+'_diff ('+var1.variables[var_name].units+') of '+str(ts_s)+'-'+str(ts_e), size=16)
            cb.set_label(var1.variables[var_name].units,size=14,rotation=0,labelpad=15)
        else:
            plt.title(var_name+'_diff ('+units+') of '+str(ts_s)+'-'+str(ts_e), size=16)
            cb.set_label(units,size=14,rotation=0,labelpad=15)
        cb.ax.tick_params(labelsize=10)

        plt.savefig('./plots/spatial_map_'+var_name+'_'+case_names[0]+'_vs_'+case_names[1]
                    + '_ts-'+str(ts_s)+'_'+str(ts_e)+'.png',dpi=1200)
        Var = None

if __name__ == "__main__":

    # ### plot AWAP Tmax
    # path      = '/g/data/w35/Shared_data/Observations/AWAP_all_variables/daily/tmax/'
    # file_name = 'AWAP_daily_tmax_1970_2019.nc'
    # file_path = path + file_name
    # var_name  = 'tmax'
    # plot_map_var(file_path, var_name)

    # ### plot LIS variables at timestep = ts
    # case_name  = 'hires_r7264'
    # file_name  = "LIS.CABLE.198212-201301.nc"
    # path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    # file_path  = path + file_name
    # # Since lon and lat in LIS contain default values, to use plt.contourf, I take lon/lat from WRF output
    # wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/wrfout_d01_2013-01-01_03:00:00"
    # d2012121   = 10524
    # ts         = d2012121+30 + 4  # 2013-01-04  #- 2012-12-01
    # var_names  = ["Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg",
    #              "Rainf_tavg","Evap_tavg","Qs_tavg","Qsb_tavg","VegT_tavg","AvgSurfT_tavg",
    #              "Albedo_inst","SoilWet_inst","ECanop_tavg","TVeg_tavg",
    #              "FWsoil_tavg","ESoil_tavg","GPP_tavg","Wind_f_inst",
    #              "Tair_f_inst", "Qair_f_inst","Psurf_f_inst","SWdown_f_inst","LWdown_f_inst"]
    # plot_map_var_ts(file_path, wrf_path, var_names, ts)



    ### plot plot_map_var_diff_period_mean
    case_names = ['free_drain_hires_r7264','hires_r7264']
    file_name  = "LIS.CABLE.198212-201301.nc"
    file_paths = []
    for case_name in case_names:
        path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
        file_path  = path + file_name
        file_paths.append(file_path)

    # Since lon and lat in LIS contain default values, to use plt.contourf, I take lon/lat from WRF output
    wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/wrfout_d01_2013-01-01_03:00:00"
    d2012121   = 10524
    ts_s       = d2012121 -365 #+ 30 + 1
    ts_e       = d2012121 # + 30 + 14  # 2013-01-04  #- 2012-12-01
    var_names  = ["Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg",
                 "Rainf_tavg","Evap_tavg","Qs_tavg","Qsb_tavg","VegT_tavg","AvgSurfT_tavg",
                 "Albedo_inst","SoilWet_inst","ECanop_tavg","TVeg_tavg",
                 "FWsoil_tavg","ESoil_tavg","GPP_tavg","Wind_f_inst",
                 "Tair_f_inst", "Qair_f_inst","Psurf_f_inst","SWdown_f_inst","LWdown_f_inst"]
    plot_map_var_diff_period_mean(file_paths, wrf_path, case_names, var_names, ts_s, ts_e)
