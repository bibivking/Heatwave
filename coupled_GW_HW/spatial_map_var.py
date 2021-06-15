#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

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

def plot_map_var_time(file_path, var_names, ts):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var = Dataset(file_path, mode='r')

    # time = var.variables['time']
    lons = var.variables['lon'] # or ['lon']
    lats = var.variables['lat']  # or ['lat']
    _FillValue = -9999.
    print(_FillValue)
    lon, lat = np.meshgrid(lons, lats)

    lon[lon == _FillValue] = np.nan
    lat[lat == _FillValue] = np.nan
    print(np.nanmean(lon))
    print(np.nanmean(lat))
    
    for var_name in var_names:

        print(var_name)
        Var = var.variables[var_name][ts]
        Var[Var == _FillValue] = np.nan

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

        # print(min(Var))
        # print(max(Var))
        # print((max(Var)-min(Var))/20.)

        # clevs = np.arange(np.min(Var),np.max(Var)+(np.max(Var)-np.max(Var))/20.,(np.max(Var)-np.max(Var))/20.)
        if len(np.shape(Var)) == 2:
            plt.contourf(lon, lat, Var[:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.jet) # , T2M_daily_avg 
            print("hi hi")
        if len(np.shape(Var)) == 3:
            plt.contourf(lon, lat, Var[0,:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.jet) # clevs, T2M_daily_avg

        # plt.title(var_name+' ('+Var.units+') ts='+str(ts), size=16)
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        # cb.set_label(Var.units,size=14,rotation=0,labelpad=15)
        cb.ax.tick_params(labelsize=10)
        # Overlay wind vectors
        # qv = plt.quiver(lon, lat, U10M[i,:,:], V10M[i,:,:], scale=420, color='k')
        plt.savefig('./plots/'+var_name+'_ts-'+str(ts)+'.png',dpi=1200)
        Var = None

if __name__ == "__main__":

    # ### plot AWAP Tmax
    # path      = '/g/data/w35/Shared_data/Observations/AWAP_all_variables/daily/tmax/'
    # file_name = 'AWAP_daily_tmax_1970_2019.nc'
    # file_path = path + file_name
    # var_name  = 'tmax'
    # plot_map_var(file_path, var_name)

    ### plot LIS variables at timestep = ts
    case_name  = 'hires_r7264'
    file_name  = "LIS.CABLE.198212-201301.nc"
    path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    file_path  = path + file_name
    tss        = [10524 + 30 + 4, 10524 + 30 + 5, 10524 + 30 + 6, 10524 + 30 + 7, 10524 + 30 + 8]  # 2013-01-04  #10524- 2012-12-01
    var_names  = ["Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg",
                 "Rainf_tavg","Evap_tavg","Qs_tavg","Qsb_tavg","VegT_tavg","AvgSurfT_tavg",
                 "Albedo_inst","SoilWet_inst","ECanop_tavg","TVeg_tavg",
                 "FWsoil_tavg","ESoil_tavg","GPP_tavg","Wind_f_inst",
                 "Tair_f_inst", "Qair_f_inst","Psurf_f_inst","SWdown_f_inst","LWdown_f_inst"]

    for ts in tss : 
        plot_map_var_time(file_path, var_names, ts)


    # case_names = ['free_drain_hires_r7264','hires_r7264']
    # file_name  = "LIS.CABLE.198212-201301.nc"
    # file_paths = []
    # for case_name in case_names:
    #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    #     file_path  = path + file_name
    #     file_paths.append(file_path)

    # file_type = 'all_in_one'
    # layer     = 6
    # var_dim   = 3
    # var_names = ["Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg","Snowf_tavg",
    #             "Rainf_tavg","Evap_tavg","Qs_tavg","Qsb_tavg","VegT_tavg","AvgSurfT_tavg",
    #             "Albedo_inst","SWE_inst","SnowDepth_inst","SoilWet_inst","ECanop_tavg","TVeg_tavg",
    #             "FWsoil_tavg","ESoil_tavg","CanopInt_inst","SnowCover_inst","GPP_tavg","Wind_f_inst",
    #             "Rainf_f_inst","Tair_f_inst", "Qair_f_inst","Psurf_f_inst","SWdown_f_inst","LWdown_f_inst"]
