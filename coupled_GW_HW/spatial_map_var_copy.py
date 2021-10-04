#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_range_diff
from common_utils import get_reverse_colormap

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

def plot_map_var_ts(is_diff, file_paths, wrf_path, case_names, var_names, tss, layer=None):

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

            if get_reverse_colormap(var_name) == None:
                cmap = plt.cm.seismic
            elif get_reverse_colormap(var_name) == True:
                cmap = plt.cm.seismic_r
            else:
                cmap = plt.cm.seismic

            if ranges == None:
                if len(np.shape(Var)) == 2:
                    plt.contourf(lon, lat, Var[:,:], transform=ccrs.PlateCarree(),cmap=cmap)
                elif len(np.shape(Var)) == 3:
                    plt.contourf(lon, lat, Var[layer,:,:], transform=ccrs.PlateCarree(),cmap=cmap)
            else:
                clevs = np.linspace(ranges[0],ranges[1], num=21)
                if len(np.shape(Var)) == 2:
                    plt.contourf(lon, lat, Var[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
                if len(np.shape(Var)) == 3:
                    plt.contourf(lon, lat, Var[layer,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)

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

            if is_diff:
                if len(np.shape(Var)) == 3: # 4-D var
                    plt.savefig('./plots/Oct2021/land_var/spatial_map_'+var_name+'_lvl-'+str(layer)+'_'+case_names[0]+'_vs_'+case_names[1]
                                + '_ts-'+str(ts)+'.png',dpi=300)
                else:
                    plt.savefig('./plots/Oct2021/land_var/spatial_map_'+var_name+'_'+case_names[0]+'_vs_'+case_names[1]
                                + '_ts-'+str(ts)+'.png',dpi=300)
            else:
                if len(np.shape(Var)) == 3:
                    plt.savefig('./plots/Oct2021/land_var/spatial_map_'+var_name+'_lvl-'+str(layer)+'_'+case_names[0]+'_ts-'+str(ts)+'.png',dpi=300)
                else:
                    plt.savefig('./plots/Oct2021/land_var/spatial_map_'+var_name+'_'+case_names[0]+'_ts-'+str(ts)+'.png',dpi=300)
            Var = None

def plot_map_var_period_mean(is_diff, file_paths, wrf_path, case_names, var_names, ts_s, ts_e,layer=None):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:

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


        if get_reverse_colormap(var_name) == None:
            cmap = plt.cm.seismic
        elif get_reverse_colormap(var_name) == True:
            cmap = plt.cm.seismic_r
        else:
            cmap = plt.cm.seismic

        if is_diff:
            if ranges == None:
                if len(np.shape(Var)) == 2:
                    plt.contourf(lon, lat, Var[:,:], transform=ccrs.PlateCarree(),cmap=cmap)
                elif len(np.shape(Var)) == 3:
                    plt.contourf(lon, lat, Var[layer,:,:], transform=ccrs.PlateCarree(),cmap=cmap)
            else:
                clevs = np.linspace(ranges[0],ranges[1], num=21)
                if len(np.shape(Var)) == 2:
                    plt.contourf(lon, lat, Var[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
                if len(np.shape(Var)) == 3:
                    plt.contourf(lon, lat, Var[layer,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
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
        if is_diff:
            if len(np.shape(Var)) == 3: # 4-D var
                plt.savefig('./plots/Oct2021/land_var/spatial_map_'+var_name+'_lvl-'+str(layer)+'_'+case_names[0]+'_vs_'+case_names[1]
                            + '_ts-'+str(ts_s)+'_'+str(ts_e)+'.png',dpi=300)
            else:
                plt.savefig('./plots/Oct2021/land_var/spatial_map_'+var_name+'_'+case_names[0]+'_vs_'+case_names[1]
                            + '_ts-'+str(ts_s)+'_'+str(ts_e)+'.png',dpi=300)
        else:
            if len(np.shape(Var)) == 3:
                plt.savefig('./plots/Oct2021/land_var/spatial_map_'+var_name+'_lvl-'+str(layer)+'_'+case_names[0]+'_ts-'+str(ts_s)+'_'+str(ts_e)+'.png',dpi=300)
            else:
                plt.savefig('./plots/Oct2021/land_var/spatial_map_'+var_name+'_'+case_names[0]+'_ts-'+str(ts_s)+'_'+str(ts_e)+'.png',dpi=300)
        Var = None

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

    # =============================== Operation ================================
    case_names = [
                    ["fd_feb2011","gw_feb2011"],
                    ["fd_jan2013","gw_jan2013"],
                    ["fd_jan2014","gw_jan2014"],
                    ["fd_feb2017","gw_feb2017"],
                    ["fd_jan2019","gw_jan2019"] ] # ["fd_feb2009","gw_feb2009"],
    file_names = [
                    "LIS.CABLE.201102-201102.d01.nc",
                    "LIS.CABLE.201301-201301.d01.nc",
                    "LIS.CABLE.201401-201401.d01.nc",
                    "LIS.CABLE.201702-201702.d01.nc",
                    "LIS.CABLE.201901-201901.d01.nc" ] # "LIS.CABLE.200902-200902.d01.nc",

    case_sum   = len(file_names)
    wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/fd_feb2009/WRF_output/wrfout_d01_2009-02-01_00:00:00"

    # ############################
    #   plot plot_map_var_ts    #
    # ############################
    # Since lon and lat in LIS contain default values, to use plt.contourf, I take lon/lat from WRF output
    is_diff    = True #False

    for case_num in np.arange(case_sum):

        file_paths = []

        for case_name in case_names[case_num]:
            print(case_name)
            path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
            file_path  = path + file_names[case_num]
            file_paths.append(file_path)

        tss        = []

        # day 1 to day 28
        for ts in np.arange(224):
            tss.append(ts)

        var_names  = ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg",
                      "Qh_tavg","Qle_tavg","Qg_tavg","Qs_tavg","Qsb_tavg",
                      "FWsoil_tavg","AvgSurfT_tavg","VegT_tavg","Tair_f_inst",
                      "Rainf_tavg", "Qair_f_inst"]

        plot_map_var_ts(is_diff, file_paths, wrf_path, case_names[case_num], var_names, tss)

        var_names  = ["SoilMoist_inst","SoilTemp_inst"]
        for layer in np.arange(0,6):
            plot_map_var_ts(is_diff, file_paths, wrf_path, case_names[case_num], var_names, tss, layer=layer)

    # ######################################
    # #   plot plot_map_var_period_mean    #
    # ######################################
    # is_diff    = True #False
    # layer      = None
    # var_names  = ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg",
    #               "Qh_tavg","Qle_tavg","Qg_tavg","Qs_tavg","Qsb_tavg",
    #               "FWsoil_tavg","AvgSurfT_tavg","VegT_tavg","Tair_f_inst",
    #               "Rainf_tavg", "Qair_f_inst"]
    #
    # for case_num in np.arange(case_sum):
    #     file_paths = []
    #     for case_name in case_names[case_num]:
    #         path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    #         file_path  = path + file_names[case_num]
    #         file_paths.append(file_path)
    #
    #     # from day 6 to day 28
    #     ts_s       = 5*8
    #     ts_e       = 28*8
    #
    #     #for layer in np.arange(0,6):
    #     plot_map_var_period_mean(is_diff, file_paths, wrf_path, case_names[case_num], var_names, ts_s, ts_e, layer=layer)

    # ==================== Old Operation (before Oct 2021) =====================
    ######################
    #   Plot AWAP Tmax   #
    ######################

    # path      = '/g/data/w35/Shared_data/Observations/AWAP_all_variables/daily/tmax/'
    # file_name = 'AWAP_daily_tmax_1970_2019.nc'
    # file_path = path + file_name
    # var_name  = 'tmax'
    # plot_map_var(file_path, var_name)

    ############################################
    #   plot LIS variables at timestep = ts    #
    ############################################
    # file_name  = "LIS.CABLE.198212-201301.nc"
    # d2012121   = 10524
    # ts         = d2012121 # d2012121+30 + 4  # 2013-01-04  #- 2012-12-01
    # # var_names  = ["Albedo_inst","Landcover_inst","Landmask_inst","Soiltype_inst","SandFrac_inst","ClayFrac_inst","SiltFrac_inst",
    # #               "SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst","Hyds_inst","Bch_inst","Sucs_inst","Elevation_inst","LAI_inst"]"GPP_tavg","VegT_tavg",

    # case_name  = 'ctl_11Jul'
    # path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    # file_path  = path + file_name
    # # Since lon and lat in LIS contain default values, to use plt.contourf, I take lon/lat from WRF output
    # wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/wrfout_d01_2013-01-01_03:00:00"
    # plot_map_var_ts(case_name,file_path, wrf_path, var_names, ts)

    # case_name  = 'free_drain_11Jul'
    # path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    # file_path  = path + file_name
    # # Since lon and lat in LIS contain default values, to use plt.contourf, I take lon/lat from WRF output
    # wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/wrfout_d01_2013-01-01_03:00:00"
    # plot_map_var_ts(case_name,file_path, wrf_path, var_names, ts)

    #############################
    #   plot plot_map_var_ts    #
    #############################

    # case_names = ['free_drain_25Jul','ctl_25Jul']
    # file_name  = "LIS.CABLE.198212-201301.nc"
    # file_paths = []
    # is_diff    = True #False
    # for case_name in case_names:
    #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    #     file_path  = path + file_name
    #     file_paths.append(file_path)
    #
    # # Since lon and lat in LIS contain default values, to use plt.contourf, I take lon/lat from WRF output
    # wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/wrfout_d01_2013-01-01_03:00:00"
    #
    # d2012121   = 10524
    # tss        = []
    #
    # # 2012-11-22 ~ 2012-12-10
    # for ts in np.arange(d2012121-10,d2012121+10):
    #     tss.append(ts)
    #
    # # 2012-12-26 ~ 2013-01-14
    # for ts in np.arange(d2012121+31-5,d2012121+31+14):
    #     tss.append(ts)
    #
    # var_names  = ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","FWsoil_tavg","Qh_tavg"]
    #             #  ["Rainf_tavg","Tair_f_inst", "Qair_f_inst"]
    # plot_map_var_ts(is_diff, file_paths, wrf_path, case_names, var_names, tss)
    #
    # var_names  = ["SoilMoist_inst","SoilTemp_inst"]
    # for layer in np.arange(0,6):
    #     plot_map_var_ts(is_diff, file_paths, wrf_path, case_names, var_names, tss, layer=layer)


    # ######################################
    # #   plot plot_map_var_period_mean    #
    # ######################################
    # case_names = ['free_drain_14Aug','ctl_14Aug'] #'free_drain_25Jul',
    # file_name  = "LIS.CABLE.198212-201301.nc"
    # file_paths = []
    # is_diff    = True #False
    # layer      = None
    # for case_name in case_names:
    #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    #     file_path  = path + file_name
    #     file_paths.append(file_path)

    # # Since lon and lat in LIS contain default values, to use plt.contourf, I take lon/lat from WRF output
    # wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/wrfout_d01_2013-01-01_03:00:00"
    # d2012121   = 10524
    # ts_s       = d2012121 -366 #+ 30 + 1
    # ts_e       = d2012121  # + 30 + 14  # 2013-01-04  #- 2012-12-01
    # var_names  = var_3D_basic_names #["WaterTableD_tavg"]#["SoilMoist_inst","SoilTemp_inst"] #var_3D_basic_names # ["SoilMoist_inst","SoilTemp_inst"]

    # #for layer in np.arange(0,6):
    # plot_map_var_period_mean(is_diff, file_paths, wrf_path, case_names, var_names, ts_s, ts_e, layer=layer)

    # ######################################
    # #   plot plot_map_var_period_mean    #
    # ######################################
    #
    # ############# duplicate ##############
    # # old run vs new run
    # case_names = ['ctl_23Aug','ctl_8Sep'] #'free_drain_25Jul',
    # file_name  = "LIS.CABLE.200001-200012.nc"#"LIS.CABLE.198212-201301.nc"
    # file_paths = []
    # is_diff    = True #False
    # layer      = None
    # for case_name in case_names:
    #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    #     file_path  = path + file_name
    #     file_paths.append(file_path)
    #
    # # Since lon and lat in LIS contain default values, to use plt.contourf, I take lon/lat from WRF output
    # wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/wrfout_d01_2013-01-01_03:00:00"
    # d2012121   = 10524
    # ts_s       = 0 #+ 30 + 1 # d2012121
    # ts_e       = 366  # + 30 + 14  # 2013-01-04  #- 2012-12-01 # d2012121
    # var_names  = var_3D_basic_names #["WaterTableD_tavg"]#["SoilMoist_inst","SoilTemp_inst"] #var_3D_basic_names # ["SoilMoist_inst","SoilTemp_inst"]
    #
    # #for layer in np.arange(0,6):
    # plot_map_var_period_mean(is_diff, file_paths, wrf_path, case_names, var_names, ts_s, ts_e, layer=layer)
