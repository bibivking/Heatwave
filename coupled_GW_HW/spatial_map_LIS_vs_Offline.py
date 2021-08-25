#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_scale_offline
from common_utils import get_reverse_colormap

def plot_map_var_offline(file_path, var_names):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var = Dataset(file_path, mode='r')

    time = var.variables['time']
    lons = var.variables['longitude'][:,:] # or ['lon']
    lats = var.variables['latitude'][:,:]  # or ['lat']
    # lon, lat = np.meshgrid(lons, lats)
    print(lons)
    print(lats)

    for var_name in var_names:
        print(var_name)
        scale, units = get_land_var_scale_offline(var_name)
        
        Var   = np.mean(var.variables[var_name],axis=0)*scale
        print(Var.shape)

        fig = plt.figure(figsize=(7,5))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([140,154,-40,-28])

        ax.coastlines(resolution="50m",linewidth=1)
        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator([140,145,150])
        gl.ylocator = mticker.FixedLocator([-40,-35,-30])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black'}
        gl.ylabel_style = {'size':10, 'color':'black'}
        # Plot windspeed

        if var_name == 'WatTable':
            clevs = np.arange(0,15,1)
        else:
            clevs = np.arange(0,3.,0.1)

        plt.contourf(lons, lats, Var, clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG)#seismic) # T2M_daily_avg 
        plt.title(var_name, size=16)
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.set_label(units,size=14,rotation=0,labelpad=15)
        cb.ax.tick_params(labelsize=10)
        plt.savefig('./plots/lis_vs_offline/spatial_map_offline_'+var_name+'.png',dpi=300)

def plot_map_var_lis(file_path, wrf_path, case_name, var_names):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:

    var = Dataset(file_path, mode='r')
   
    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    for var_name in var_names:

        print(var_name)
        scale, units = get_land_var_scale(var_name)

        if scale == -273.15:
            Var   = np.mean(var.variables[var_name],axis=0)
        else:
            Var   = np.mean(var.variables[var_name],axis=0)*scale
        print(Var.shape)

        # Make plots
        fig = plt.figure(figsize=(7,5))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([140,154,-40,-28])
        ax.coastlines(resolution="50m",linewidth=1)

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator([140,145,150])
        gl.ylocator = mticker.FixedLocator([-40,-35,-30])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black'}
        gl.ylabel_style = {'size':10, 'color':'black'}

        # if get_reverse_colormap(var_name) == None:
        #     cmap = plt.cm.seismic
        # elif get_reverse_colormap(var_name) == True:
        #     cmap = plt.cm.seismic_r
        # else:
        #     cmap = plt.cm.seismic
        if var_name == 'WaterTableD_tavg':
            clevs = np.arange(0,15,1)
        else:
            clevs = np.arange(0,3,0.1)

        plt.contourf(lon, lat, Var[:,:],clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG)#seismic)

        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

        if units == None:
            units_string = var.variables[var_name].units
        else:
            units_string = units

        plt.title(var_name+' ('+units_string+')', size=16)

        cb.ax.tick_params(labelsize=10)
        cb.set_label(units_string,size=14,rotation=270,labelpad=15)

        plt.savefig('./plots/lis_vs_offline/spatial_map_lis_'+var_name+'_'+case_name+'.png',dpi=300)
        Var = None

def statistic_bar(file_path_offline, file_path_lis, wrf_path):

    off_stats   = np.zeros(7)
    lis_stats   = np.zeros(7)
    scale       = 3600*24.*366.

    # Read offline
    offline = Dataset(file_path_offline, mode='r')

    def_off = offline.variables['Rainf']._FillValue
    
    rain_off_tmp = np.mean(offline.variables['Rainf'],axis=0)*scale  
    off_stats[0] = np.nanmean(np.where(rain_off_tmp == def_off, np.nan, rain_off_tmp))
    
    evap_off_tmp = np.mean(offline.variables['Evap'],axis=0)*scale  
    off_stats[1] = np.nanmean(np.where(evap_off_tmp == def_off, np.nan, evap_off_tmp))

    tveg_off_tmp = np.mean(offline.variables['TVeg'],axis=0)*scale  
    off_stats[2] = np.nanmean(np.where(tveg_off_tmp == def_off, np.nan, tveg_off_tmp))
    
    esoil_off_tmp = np.mean(offline.variables['ESoil'],axis=0)*scale  
    off_stats[3]  = np.nanmean(np.where(esoil_off_tmp == def_off, np.nan, esoil_off_tmp))

    ecanop_off_tmp = np.mean(offline.variables['ECanop'],axis=0)*scale  
    off_stats[4]   = np.nanmean(np.where(ecanop_off_tmp == def_off, np.nan, ecanop_off_tmp))

    qs_off_tmp     = np.mean(offline.variables['Qs'],axis=0)*scale  
    off_stats[5]   = np.nanmean(np.where(qs_off_tmp == def_off, np.nan, qs_off_tmp))

    qsb_off_tmp    = np.mean(offline.variables['Qsb'],axis=0)*scale  
    off_stats[6]   = np.nanmean(np.where(qsb_off_tmp == def_off, np.nan, qsb_off_tmp))    

    # Read lis-cable
    lis     = Dataset(file_path_lis, mode='r')
    wrf     = Dataset(wrf_path,  mode='r')

    # mask: keep SE Aus
    lon     = wrf.variables['XLONG'][0,:,:]
    lat     = wrf.variables['XLAT'][0,:,:]
    mask    = (lon >= 140) & (lon <= 154) & (lat >= -40) & (lat <= -28) 

    def_lis = lis.variables['Rainf_tavg']._FillValue
    
    rain_lis_tmp = np.mean(lis.variables['Rainf_tavg'],axis=0)*scale
    lis_stats[0] = np.nanmean(np.where(rain_lis_tmp[mask] == def_lis, np.nan, rain_lis_tmp[mask]))
    
    evap_lis_tmp = np.mean(lis.variables['Evap_tavg'],axis=0)*scale
    lis_stats[1] = np.nanmean(np.where(evap_lis_tmp[mask] == def_lis, np.nan, evap_lis_tmp[mask]))

    tveg_lis_tmp = np.mean(lis.variables['TVeg_tavg'],axis=0)*scale  
    lis_stats[2] = np.nanmean(np.where(tveg_lis_tmp[mask] == def_lis, np.nan, tveg_lis_tmp[mask]))
    
    esoil_lis_tmp = np.mean(lis.variables['ESoil_tavg'],axis=0)*scale  
    lis_stats[3]  = np.nanmean(np.where(esoil_lis_tmp[mask] == def_lis, np.nan, esoil_lis_tmp[mask]))

    ecanop_lis_tmp = np.mean(lis.variables['ECanop_tavg'],axis=0)*scale  
    lis_stats[4]   = np.nanmean(np.where(ecanop_lis_tmp[mask] == def_lis, np.nan, ecanop_lis_tmp[mask]))

    qs_lis_tmp     = np.mean(lis.variables['Qs_tavg'],axis=0)*scale  
    lis_stats[5]   = np.nanmean(np.where(qs_lis_tmp[mask] == def_lis, np.nan, qs_lis_tmp[mask]))

    qsb_lis_tmp    = np.mean(lis.variables['Qsb_tavg'],axis=0)*scale  
    lis_stats[6]   = np.nanmean(np.where(qsb_lis_tmp[mask] == def_lis, np.nan, qsb_lis_tmp[mask]))    

    # plotting
    labels = ['Rain', 'Evap', 'TVeg', 'Esoil', 'Ecan','Qs','Qsb']

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, off_stats, width, label='Off')
    rects2 = ax.bar(x + width/2, lis_stats, width, label='LIS')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('mm/yr')
    ax.set_title('Water Balance')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)

    fig.tight_layout()

    plt.savefig('./plots/lis_vs_offline/water_balance_lis_vs_off_statistic.png',dpi=300)


if __name__ == "__main__":

    var_LIS_names      = ['Rainf_tavg','Evap_tavg','ESoil_tavg','ECanop_tavg','TVeg_tavg',
                          'Qs_tavg','Qsb_tavg','WaterTableD_tavg','Qle_tavg','Qh_tavg','Qg_tavg']

    var_offline_names  = ['Rainf','Evap','ESoil','ECanop','TVeg','Qs','Qsb','WatTable','Qle','Qh','Qg']

    ##########################
    #   Plot offline CABLE   #
    ##########################

    path      = '/g/data/w35/mm3972/model/cable/runs/test_para_chg_dpt/uniform_6layer/outputs/'
    file_name = 'cable_out_2000_SE_Aus.nc'
    file_path = path + file_name
    plot_map_var_offline(file_path, var_offline_names)

    ##############################
    #   plot plot_map_var_lis    # 
    ##############################
    case_name  = 'ctl_23Aug'
    file_name  = "LIS.CABLE.200001-200012.nc"
    layer      = None
    path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_23Aug/LIS_output/"
    file_path  = path + file_name
    wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_14Aug/WRF_output_copy/wrfout_d01_2013-01-01_03:00:00"
    var_names  = var_LIS_names 

    plot_map_var_lis(file_path, wrf_path, case_name, var_names)

    ##############################
    #   plot statistic_bar    # 
    ##############################
    file_path_offline = '/g/data/w35/mm3972/model/cable/runs/test_para_chg_dpt/uniform_6layer/outputs/cable_out_2000_SE_Aus.nc'
    file_path_lis     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_23Aug/LIS_output/LIS.CABLE.200001-200012.nc'
    wrf_path          = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_14Aug/WRF_output_copy/wrfout_d01_2013-01-01_03:00:00"
    statistic_bar(file_path_offline, file_path_lis, wrf_path)