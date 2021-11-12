#!/usr/bin/python

'''
Plot spitial map of land diagnosis and parameters from LIS-CABLE
1. per time step
2. time period average
'''

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_range_diff
from common_utils import *

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
        # plt.contourf(lon, lat, ws_daily_avg, clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.seismic)
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
        plt.contourf(lon, lat, Var[i,:,:], clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.seismic) # T2M_daily_avg
        plt.title('AWAP Tmax ' + str(i-15675), size=16)
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.set_label('C',size=14,rotation=0,labelpad=15)
        cb.ax.tick_params(labelsize=10)
        # Overlay wind vectors
        # qv = plt.quiver(lon, lat, U10M[i,:,:], V10M[i,:,:], scale=420, color='k')
        plt.savefig('./plots/'+var_name+'-'+str(i-15675)+'.png',dpi=300)

def plot_map_var_ts(is_diff, file_paths, wrf_path, var_names, tss, layer=None, message=None):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var1 = Dataset(file_paths[0], mode='r')
    if is_diff:
        var2 = Dataset(file_paths[1], mode='r')

    wrf = Dataset(wrf_path,  mode='r')

    # use WRF output's lat & lon, since LIS output has default value
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]
    for ts in tss:
        for var_name in var_names:
            print(var_name)
            scale, units = get_land_var_scale(var_name)
            if is_diff:
                ranges = get_land_var_range_diff(var_name)
                if scale == -273.15:
                    Var = var2.variables[var_name][ts] - var1.variables[var_name][ts]
                else:
                    Var = (var2.variables[var_name][ts] - var1.variables[var_name][ts])*scale
            else:
                if scale == -273.15:
                    Var = var1.variables[var_name][ts] + scale
                else:
                    Var = var1.variables[var_name][ts]*scale

            # Make plots
            fig = plt.figure(figsize=(9,5))
            ax = plt.axes(projection=ccrs.PlateCarree())
            # ax.set_extent([110,155,-45,-10])
            ax.set_extent([135,155,-40,-25])
            ax.coastlines(resolution="50m",linewidth=1)

            # Add gridlines
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
            gl.xlabels_top  = False
            gl.ylabels_right= False
            gl.xlines       = True
            # gl.xlocator     = mticker.FixedLocator([110,115,120,125,130,135,140,145,150,155])
            # gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15,-10])
            gl.xlocator     = mticker.FixedLocator([135,140,145,150,155])
            gl.ylocator     = mticker.FixedLocator([-40,-35,-30,-25])
            gl.xformatter   = LONGITUDE_FORMATTER
            gl.yformatter   = LATITUDE_FORMATTER
            gl.xlabel_style = {'size':10, 'color':'black'}
            gl.ylabel_style = {'size':10, 'color':'black'}

            if get_reverse_colormap(var_name) == None:
                cmap = plt.cm.seismic
            elif get_reverse_colormap(var_name) == True:
                cmap = plt.cm.seismic_r
            else:
                cmap = plt.cm.seismic

            if is_diff:
                if var_name in ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg","Rainf_tavg"]:
                    clevs = np.linspace(-5.,5., num=11)
                elif var_name in ["Qh_tavg","Qle_tavg","Qg_tavg"]:
                    clevs = np.linspace(-80.,80., num=11)
                elif var_name in ["SoilMoist_inst","Qair_f_inst"]:
                    clevs = np.linspace(-0.3,0.3, num=11)
                elif var_name in ["AvgSurfT_tavg","VegT_tavg","Tair_f_inst","SoilTemp_inst"]:
                    clevs = np.linspace(-4.,4., num=11)
                elif var_name in ["FWsoil_tavg"]:
                    clevs = np.linspace(-0.5,0.5, num=11)
                elif var_name in ["WaterTableD_tavg"]:
                    clevs = np.linspace(-10.,10., num=21)
                if len(np.shape(Var)) == 2:
                    plt.contourf(lon, lat, Var[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap, extend='both')
                if len(np.shape(Var)) == 3:
                    plt.contourf(lon, lat, Var[layer,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap, extend='both')
            else:
                if len(np.shape(Var)) == 2:
                    plt.contourf(lon, lat, Var[:,:], transform=ccrs.PlateCarree(),cmap=cmap, extend='both')
                if len(np.shape(Var)) == 3:
                    plt.contourf(lon, lat, Var[layer,:,:], transform=ccrs.PlateCarree(),cmap=cmap, extend='both')

            cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

            if units == None:
                units_string = var1.variables[var_name].units
            else:
                units_string = units

            if is_diff:
                plt.title(var_name+'_diff ('+units_string+') ts= '+str(ts), size=16)
            else:
                plt.title(var_name+' ('+units_string+') ts= '+str(ts), size=16)

            cb.ax.tick_params(labelsize=10)
            cb.set_label(units_string,size=14,rotation=270,labelpad=15)

            if message == None:
                txt = var_name
            else:
                txt = message + "_" + var_name

            if len(np.shape(Var)) == 3:
                plt.savefig('./plots/19Oct/land_var/spatial_map_'+txt+'_lvl-'+str(layer)+'_ts-'+str(ts)+'.png',dpi=300)
            else:
                plt.savefig('./plots/19Oct/land_var/spatial_map_'+txt+'_ts-'+str(ts)+'.png',dpi=300)
            Var = None

def plot_map_var_period_mean(is_diff, file_paths, wrf_path, var_names, ts_s, ts_e, layer=None, message=None):

    var1 = Dataset(file_paths[0], mode='r')
    if is_diff:
        var2 = Dataset(file_paths[1], mode='r')

    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    for var_name in var_names:

        print(var_name)
        scale, units = get_land_var_scale(var_name)
        if is_diff:
            ranges = get_land_var_range_diff(var_name)
            if scale == -273.15:
                Var   = ( np.mean(var2.variables[var_name][ts_s:ts_e],axis=0)
                        - np.mean(var1.variables[var_name][ts_s:ts_e],axis=0))
            else:
                Var   = ( np.mean(var2.variables[var_name][ts_s:ts_e],axis=0)
                        - np.mean(var1.variables[var_name][ts_s:ts_e],axis=0))*scale
        else:
            if scale == -273.15:
                Var   = np.mean(var1.variables[var_name][ts_s:ts_e],axis=0)
            else:
                Var   = np.mean(var1.variables[var_name][ts_s:ts_e],axis=0)*scale
        # Make plots
        fig = plt.figure(figsize=(9,5))
        ax = plt.axes(projection=ccrs.PlateCarree())
        # ax.set_extent([110,155,-45,-10])
        ax.set_extent([135,155,-40,-25])
        ax.coastlines(resolution="50m",linewidth=1)

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = True
        # gl.xlocator     = mticker.FixedLocator([110,115,120,125,130,135,140,145,150,155])
        # gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15,-10])
        gl.xlocator     = mticker.FixedLocator([135,140,145,150,155])
        gl.ylocator     = mticker.FixedLocator([-40,-35,-30,-25])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black'}
        gl.ylabel_style = {'size':10, 'color':'black'}

        if get_reverse_colormap(var_name) == None:
            cmap = plt.cm.seismic
        elif get_reverse_colormap(var_name) == True:
            cmap = plt.cm.seismic_r
        else:
            cmap = plt.cm.seismic

        if is_diff:
            if var_name in ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg","Rainf_tavg"]:
                clevs = np.linspace(-5.,5., num=11)
            elif var_name in ["Qh_tavg","Qle_tavg","Qg_tavg"]:
                clevs = np.linspace(-50.,50., num=11)
            elif var_name in ["SoilMoist_inst","Qair_f_inst"]:
                clevs = np.linspace(-0.3,0.3, num=11)
            elif var_name in ["AvgSurfT_tavg","VegT_tavg","Tair_f_inst","SoilTemp_inst"]:
                clevs = np.linspace(-4.,4., num=11)
            elif var_name in ["FWsoil_tavg"]:
                clevs = np.linspace(-0.5,0.5, num=11)
            elif var_name in ["WaterTableD_tavg"]:
                clevs = np.linspace(-10.,10., num=21)
            if len(np.shape(Var)) == 2:
                plt.contourf(lon, lat, Var[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap, extend='both')
            if len(np.shape(Var)) == 3:
                plt.contourf(lon, lat, Var[layer,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap, extend='both')
        else:
            if len(np.shape(Var)) == 2:
                plt.contourf(lon, lat, Var[:,:], transform=ccrs.PlateCarree(),cmap=cmap)
            elif len(np.shape(Var)) == 3:
                plt.contourf(lon, lat, Var[layer,:,:], transform=ccrs.PlateCarree(),cmap=cmap)

        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

        if units == None:
            units_string = var1.variables[var_name].units
        else:
            units_string = units

        if is_diff:
            plt.title(var_name+'_diff ('+units_string+') of '+str(ts_s)+'-'+str(ts_e), size=16)
        else:
            plt.title(var_name+' ('+units_string+') of '+str(ts_s)+'-'+str(ts_e), size=16)

        cb.ax.tick_params(labelsize=10)
        cb.set_label(units_string,size=14,rotation=270,labelpad=15)
        if message == None:
            txt = var_name
        else:
            txt = message + "_" + var_name

        if len(np.shape(Var)) == 3:
            plt.savefig('./plots/19Oct/land_var/spatial_map_'+txt+'_lvl-'+str(layer)+'_ts-'+str(ts_s)+'_'+str(ts_e)+'.png',dpi=300)
        else:
            plt.savefig('./plots/19Oct/land_var/spatial_map_'+txt+'_ts-'+str(ts_s)+'_'+str(ts_e)+'.png',dpi=300)

        Var = None

def plot_spatial_land(file_paths, wrf_path, var_name, time_s, time_e, seconds=None,loc_lat=None, loc_lon=None, message=None):

    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    file1 = Dataset(file_paths[0], mode='r')
    Time = nc.num2date(file1.variables['time'][:],file1.variables['time'].units,
            only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time = Time - datetime(2000,1,1,0,0,0)

    if var_name == "Rnet":
        Var1 = file1.variables["Qle_tavg"][:] + file1.variables["Qh_tavg"][:] + file1.variables["Qg_tavg"][:]
    elif var_name == "EF":
        Var1 = np.where( (file1.variables["Qle_tavg"][:]+file1.variables["Qh_tavg"][:]) != 0,
                file1.variables["Qle_tavg"][:]/(file1.variables["Qle_tavg"][:]+file1.variables["Qh_tavg"][:]),
                np.nan)
    else:
        Var1 = file1.variables[var_name]

    var1 = spital_var(time,Var1,time_s,time_e,seconds)

    if len(file_paths) > 1:
        file2 = Dataset(file_paths[1], mode='r')
        if var_name == "Rnet":
            Var2 = file2.variables["Qle_tavg"][:] + file2.variables["Qh_tavg"][:] + file2.variables["Qg_tavg"][:]
        elif var_name == "EF":
            Var2 = np.where( (file2.variables["Qle_tavg"][:]+file2.variables["Qh_tavg"][:]) != 0,
                    file2.variables["Qle_tavg"][:]/(file2.variables["Qle_tavg"][:]+file2.variables["Qh_tavg"][:]),
                    np.nan)
        else:
            Var2 = file2.variables[var_name]
        var2 = spital_var(time,Var2,time_s,time_e,seconds)
        var  = var2 - var1
    else:
        var  = var1

    # Make plots
    fig = plt.figure(figsize=(7,5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    # ax.set_extent([135,155,-40,-25])
    ax.coastlines(resolution="50m",linewidth=1)

    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = True
    gl.xlocator     = mticker.FixedLocator([110,115,120,125,130,135,140,145,150,155])
    gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15,-10])
    gl.xformatter   = LONGITUDE_FORMATTER
    gl.yformatter   = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':'black'}
    gl.ylabel_style = {'size':10, 'color':'black'}

    if get_reverse_colormap(var_name) == None:
        cmap = plt.cm.seismic
    elif get_reverse_colormap(var_name) == True:
        cmap = plt.cm.seismic_r
    else:
        cmap = plt.cm.seismic_r

    if len(file_paths) > 1:
        if var_name == 'EF':
            clevs = np.linspace(-0.5,0.5, num=21)
        else:
            clevs = np.linspace(-50.,50., num=21)
    else:
        clevs = np.arange(np.nanmin(var), np.nanmax(var), 21)

    plt.contourf(lon, lat, var, levels=clevs[clevs!=0], transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
    cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

    plt.title(var_name, size=16)
    cb.ax.tick_params(labelsize=10)

    if message == None:
        message = var_name
    else:
        message = message+"_"+var_name

    plt.savefig('./plots/5Nov/land_var/3Nov/spatial_map_'+message+'.png',dpi=300)

if __name__ == "__main__":

    # =============================== Variables ================================
    var_3D_names =  [ "Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg","Snowf_tavg",
                      "Rainf_tavg","Evap_tavg","Qs_tavg","Qsb_tavg","VegT_tavg","AvgSurfT_tavg",
                      "Albedo_inst","SWE_inst","SnowDepth_inst","SoilWet_inst","ECanop_tavg","TVeg_tavg",
                      "FWsoil_tavg","ESoil_tavg","CanopInt_inst","SnowCover_inst","GPP_tavg","Wind_f_inst",
                      "Rainf_f_inst","Tair_f_inst", "Qair_f_inst","Psurf_f_inst","SWdown_f_inst","LWdown_f_inst"]

    var_landinfo_3D_names =  ["Landmask_inst","Landcover_inst","Soiltype_inst","SandFrac_inst","ClayFrac_inst","SiltFrac_inst",
                              "SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst","Hyds_inst","Bch_inst","Sucs_inst",
                              "Elevation_inst","LAI_inst"]

    var_4D_names =  ["RelSMC_inst","SoilMoist_inst","SoilTemp_inst","SmLiqFrac_inst","SmFrozFrac_inst"]

    var_3D_basic_names = ['Evap_tavg',"ESoil_tavg","ECanop_tavg",'TVeg_tavg',"FWsoil_tavg","Qle_tavg","Qh_tavg","Qg_tavg","VegT_tavg","WaterTableD_tavg"]

    var_energy_names = ["Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg","Rnet","EF"]

    # =============================== Operation ================================
    case_name  = "hw2009_3Nov"

    wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/ensemble_avg/wrfout_20090122-20090213_gw"

    cpl_land_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'
    cpl_land_file_gw  = cpl_land_file + '/LIS.CABLE.20090122-20090213_gw.nc'  # land output of wrf-cable run
    cpl_land_file_fd  = cpl_land_file + '/LIS.CABLE.20090122-20090213_fd.nc'  # land output of wrf-cable run

    file_paths        = [cpl_land_file_fd,cpl_land_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw
    seconds           = [8.*60.*60.,20.*60.*60.]
    for var_name in var_energy_names:

        time_s = datetime(2009,1,28,14,0,0,0)
        time_e = datetime(2009,2,9,13,59,0,0)

        if len(file_paths) > 1:
            message = "Land_nighttime_GW-FD_"+str(time_s)
        else:
            message = "Land_nighttime_GW_"+str(time_s)

        plot_spatial_land(file_paths, wrf_path, var_name, time_s, time_e, seconds=seconds, message=message)

    # ============================= Old Decks ===============================
    # # ############################
    # #   plot plot_map_var_ts    #
    # # ############################
    # # Since lon and lat in LIS contain default values, to use plt.contourf, I take lon/lat from WRF output
    # is_diff    = False #False
    #
    #     tss        = [0]
    #     message    = case_names[case_num]+"_GW"
    #
    #     var_names  = ["Evap_tavg","TVeg_tavg","ESoil_tavg",
    #                     "Qh_tavg","Qle_tavg","FWsoil_tavg","AvgSurfT_tavg","VegT_tavg","Tair_f_inst",
    #                     "Rainf_tavg", "Qair_f_inst","WaterTableD_tavg"]
    #
    #     plot_map_var_ts(is_diff, file_paths, wrf_path, var_names, tss, message=message)
    #
    #     var_names  = ["SoilMoist_inst","SoilTemp_inst"]
    #     for layer in np.arange(0,6):
    #         plot_map_var_ts(is_diff, file_paths, wrf_path, var_names, tss, layer=layer, message=message)
    #
    # ######################################
    # #   plot plot_map_var_period_mean    #
    # ######################################
    #
    # is_diff    = False #False
    # layer      = None
    # var_names  = ["Evap_tavg","TVeg_tavg","ESoil_tavg",
    #               "Qh_tavg","Qle_tavg","FWsoil_tavg","AvgSurfT_tavg","VegT_tavg","Tair_f_inst",
    #               "Rainf_tavg", "Qair_f_inst","WaterTableD_tavg"]
    # # var_names  = ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg",
    # #               "Qh_tavg","Qle_tavg","Qg_tavg","Qs_tavg","Qsb_tavg",
    # #               "FWsoil_tavg","AvgSurfT_tavg","VegT_tavg","Tair_f_inst",
    # #               "Rainf_tavg", "Qair_f_inst","WaterTableD_tavg"]
    #
    # ts_s       = [ 6*48, 6*48, 6*48, 6*48]
    # ts_e       = [ 17*48, 13*48, 14*48, 20*48 ]
    #
    # for case_num in np.arange(case_sum):
    #     print(ts_s)
    #     print(ts_e)
    #     file_paths = []
    #
    #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_names[case_num]+"/ensemble_avg/"
    #     file_path  = path + file_names[case_num]+"_gw.nc"
    #     file_paths.append(file_path)
    #     file_path  = path + file_names[case_num]+"_gw.nc"
    #     file_paths.append(file_path)
    #     print(file_paths)
    #
    #     message    = case_names[case_num]+"_GW"
    #     plot_map_var_period_mean(is_diff, file_paths, wrf_path, var_names, ts_s[case_num], ts_e[case_num], layer=layer, message=message)
    #
    #     var_names  = ["SoilMoist_inst","SoilTemp_inst"]
    #     for layer in np.arange(0,6):
    #         plot_map_var_period_mean(is_diff, file_paths, wrf_path, var_names, ts_s[case_num], ts_e[case_num], layer=layer, message=message)
    #
    # # ############################
    # #   plot plot_map_var_ts    #
    # # ############################
    # # Since lon and lat in LIS contain default values, to use plt.contourf, I take lon/lat from WRF output
    # is_diff    = False #False
    #
    # for case_num in np.arange(case_sum):
    #
    #     file_paths = []
    #
    #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_names[case_num]+"/ensemble_avg/"
    #     file_path  = path + file_names[case_num]+"_gw.nc"
    #     file_paths.append(file_path)
    #     file_path  = path + file_names[case_num]+"_gw.nc"
    #     file_paths.append(file_path)
    #     print(file_paths)
    #
    #     tss        = [0]
    #     message    = case_names[case_num]+"_GW"
    #
    #     var_names  = ["Evap_tavg","TVeg_tavg","ESoil_tavg",
    #                     "Qh_tavg","Qle_tavg","FWsoil_tavg","AvgSurfT_tavg","VegT_tavg","Tair_f_inst",
    #                     "Rainf_tavg", "Qair_f_inst","WaterTableD_tavg"]
    #
    #     plot_map_var_ts(is_diff, file_paths, wrf_path, var_names, tss, message=message)
    #
    #     var_names  = ["SoilMoist_inst","SoilTemp_inst"]
    #     for layer in np.arange(0,6):
    #         plot_map_var_ts(is_diff, file_paths, wrf_path, var_names, tss, layer=layer, message=message)
    #
    # ######################################
    # #   plot plot_map_var_period_mean    #
    # ######################################
    #
    # is_diff    = False #False
    # layer      = None
    # var_names  = ["Evap_tavg","TVeg_tavg","ESoil_tavg",
    #               "Qh_tavg","Qle_tavg","FWsoil_tavg","AvgSurfT_tavg","VegT_tavg","Tair_f_inst",
    #               "Rainf_tavg", "Qair_f_inst","WaterTableD_tavg"]
    # # var_names  = ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg",
    # #               "Qh_tavg","Qle_tavg","Qg_tavg","Qs_tavg","Qsb_tavg",
    # #               "FWsoil_tavg","AvgSurfT_tavg","VegT_tavg","Tair_f_inst",
    # #               "Rainf_tavg", "Qair_f_inst","WaterTableD_tavg"]
    #
    # ts_s       = [ 6*48, 6*48, 6*48, 6*48]
    # ts_e       = [ 17*48, 13*48, 14*48, 20*48 ]
    #
    # for case_num in np.arange(case_sum):
    #     print(ts_s)
    #     print(ts_e)
    #     file_paths = []
    #
    #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_names[case_num]+"/ensemble_avg/"
    #     file_path  = path + file_names[case_num]+"_gw.nc"
    #     file_paths.append(file_path)
    #     file_path  = path + file_names[case_num]+"_gw.nc"
    #     file_paths.append(file_path)
    #     print(file_paths)
    #
    #     message    = case_names[case_num]+"_GW"
    #     plot_map_var_period_mean(is_diff, file_paths, wrf_path, var_names, ts_s[case_num], ts_e[case_num], layer=layer, message=message)
    #
    #     var_names  = ["SoilMoist_inst","SoilTemp_inst"]
    #     for layer in np.arange(0,6):
    #         plot_map_var_period_mean(is_diff, file_paths, wrf_path, var_names, ts_s[case_num], ts_e[case_num], layer=layer, message=message)
