#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_scale_offline
from common_utils import get_reverse_colormap

def plot_map_var_offline(file_path, var_names, is_lnd=False, layer=None, cp_path=None):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var = Dataset(file_path, mode='r')

    time = var.variables['time']
    lon  = var.variables['lon'] # or ['lon']
    lat  = var.variables['lat'][::-1] # or ['lat']
    lons, lats = np.meshgrid(lon, lat)
    print(lons)
    print(lats)

    for var_name in var_names:
        print(var_name)
        scale, units = get_land_var_scale_offline(var_name)
        var_tmp = var.variables[var_name][:,::-1,:]
        
        if cp_path!= None:
            var = Dataset(cp_path, mode='r')
            cp_tmp = var.variables["cp"][:,::-1,:]
            var_tmp = var_tmp + cp_tmp
        var_tmp[3:2928:4,:,:] = var_tmp[3:2928:4,:,:] - var_tmp[2:2928:4,:,:]
        var_tmp[2:2928:4,:,:] = var_tmp[2:2928:4,:,:] - var_tmp[1:2928:4,:,:]
        var_tmp[1:2928:4,:,:] = var_tmp[1:2928:4,:,:] - var_tmp[0:2928:4,:,:]

        Var   = np.sum(var_tmp,axis=0)*1000./366.
        
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

        if is_lnd:
            print("is_lnd"+str(is_lnd))
            if var_name in ["iveg","isoil"]:
                clevs = np.linspace(0.5,15.5, num=16)
            elif var_name in ["Albedo","sfc","ssat","swilt"]:
                clevs = np.linspace(0.0,0.5, num=21)
            elif var_name in ["sand","clay","silt"]:
                clevs = np.linspace(0.0,1.0, num=21)
            elif var_name in ["LAI"]:
                clevs = np.linspace(0.,8., num=16)
            elif var_name in ["bch"]:
                clevs = np.linspace(0.,13., num=26)
            elif var_name in ["hyds"]:
                clevs = np.linspace(0.,0.05, num=21)
            elif var_name in ["sucs"]:
                clevs = np.linspace(0.,1000., num=21)
            else:
                clevs = np.linspace(np.min(Var),np.max(Var), num=21)
        else:
            if var_name == 'WatTable':
                clevs = np.arange(0,15,1)
            elif var_name == 'SoilMoist':
                clevs = np.arange(0,0.5,0.05)
            else:
                clevs = np.arange(0,5.,0.5)
        clevs = None
        clevs = np.arange(0,5.,0.5)
        if len(np.shape(Var)) == 2:
            plt.contourf(lons, lats, Var, clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG)# hsv) #seismic) # T2M_daily_avg #BrBG
        else:
            plt.contourf(lons, lats, Var[layer,:,:], clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG) #seismic) # T2M_daily_avg  

        plt.title(var_name, size=16)
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.set_label(units,size=14,rotation=270,labelpad=15)
        cb.ax.tick_params(labelsize=10)

        if cp_path != None:
            message = var_name+"_add_cp"
        else:
            message = var_name

        if len(np.shape(Var)) == 2:
            plt.savefig('./plots/lis_vs_offline/spatial_map_offline_'+message+'.png',dpi=300) 
        else:
            plt.savefig('./plots/lis_vs_offline/spatial_map_offline_'+message+'_lvl-'+str(layer)+'.png',dpi=300) 
            
def plot_map_var_lis(file_path, wrf_path, case_name, var_names, is_lnd=False, layer=None):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:

    var = Dataset(file_path, mode='r')
   
    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    for var_name in var_names:

        print(var_name)
        scale, units = get_land_var_scale(var_name)

        Var   = np.sum(var.variables[var_name],axis=0)*24*60*60./366.
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

        if is_lnd:
            print("is_lnd"+str(is_lnd))
            if var_name in ["Landcover_inst","Soiltype_inst"]:
                clevs = np.linspace(0.5,15.5, num=16)
            elif var_name in ["Albedo_inst","SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst"]:
                clevs = np.linspace(0.0,0.5, num=21)
            elif var_name in ["SandFrac_inst","ClayFrac_inst","SiltFrac_inst"]:
                clevs = np.linspace(0.0,1.0, num=21)
            elif var_name in ["LAI_inst"]:
                clevs = np.linspace(0.,8., num=17)
            elif var_name in ["Bch_inst"]:
                clevs = np.linspace(0.,13., num=27)
            elif var_name in ["Hyds_inst"]:
                clevs = np.linspace(0.,0.05, num=21)
            elif var_name in ["Sucs_inst"]:
                clevs = np.linspace(0.,1000., num=21)
            else:
                clevs = np.linspace(np.min(Var),np.max(Var), num=21)            
        else:        
            if var_name == 'WaterTableD_tavg':
                clevs = np.arange(0,15,1)
            elif var_name == 'SoilMoist_inst':
                clevs = np.arange(0,0.5,0.05)
            else:
                clevs = np.arange(0,5,0.5)
        clevs = None
        clevs = np.arange(0,5.,0.5)
        print(Var.shape)
        if len(np.shape(Var)) == 2:
            plt.contourf(lon, lat, Var[:,:], clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG)#seismic) 
        else:
            plt.contourf(lon, lat, Var[layer,:,:], clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG)#seismic)
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

        if units == None:
            units_string = var.variables[var_name].units
        else:
            units_string = units

        plt.title(var_name+' ('+units_string+')', size=16)

        cb.ax.tick_params(labelsize=10)
        cb.set_label(units_string,size=14,rotation=270,labelpad=15)
        if len(np.shape(Var)) == 2:
            plt.savefig('./plots/lis_vs_offline/spatial_map_lis_'+var_name+'_'+case_name+'.png',dpi=300)
        else:
            plt.savefig('./plots/lis_vs_offline/spatial_map_lis_'+var_name+'_'+case_name+'_lvl-'+str(layer)+'.png',dpi=300)
            

        Var = None


if __name__ == "__main__":

    var_LIS_names      = ['Rainf_tavg','Evap_tavg','ESoil_tavg','ECanop_tavg','TVeg_tavg',
                          'Qs_tavg','Qsb_tavg','WaterTableD_tavg','Qle_tavg','Qh_tavg','Qg_tavg']

    var_offline_names  = ['Rainf','Evap','ESoil','ECanop','TVeg','Qs','Qsb','WatTable','Qle','Qh','Qg']

    var_LIS_soil_names      = ["SoilMoist_inst"]

    var_offline_soil_names  = ["SoilMoist"]

    var_LIS_landinfo_names =  ["SandFrac_inst","ClayFrac_inst","SiltFrac_inst",
                              "SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst","Hyds_inst","Bch_inst","Sucs_inst",
                              "LAI_inst"]#,"Albedo_inst","Elevation_inst"] "Landcover_inst","Soiltype_inst",

    var_offline_landinfo_names = ["sand","clay","silt",
                                  "sfc","ssat","swilt","hyds","bch","sucs", "LAI"] #"iveg","isoil",'elev',"Albedo"
    ##########################
    #   Plot offline CABLE   #
    ##########################
    is_lnd    = True
    add_cp    = True # add convective precipiation
    path      = '/g/data/w35/mm3972/model/cable/runs/test_para_chg_dpt/uniform_6layer/outputs/'
    file_name = 'cable_out_2000_SE_Aus.nc'
    file_path = '/g/data/w35/mm3972/scripts/Heatwave/coupled_GW_HW/nc_file/tp_3hrs_ERAI_historical_fc-sfc_200001_200012.nc'
    if add_cp:
        cp_path = '/g/data/w35/mm3972/scripts/Heatwave/coupled_GW_HW/nc_file/cp_3hrs_ERAI_historical_fc-sfc_200001_200012.nc'

    var_names = ['tp'] #var_offline_soil_names # var_offline_names # var_offline_soil_names
    # for layer in np.arange(6):
    plot_map_var_offline(file_path, var_names, is_lnd=is_lnd, cp_path=cp_path)#, layer=layer)

    ##############################
    #   plot plot_map_var_lis    # 
    ##############################    
    is_lnd     = True
    case_name  = 'ctl_8Sep'
    file_name  = "LIS.CABLE.200001-200012.nc"
    layer      = None
    path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_23Aug/LIS_output_2000-2004/"
    file_path  = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_8Sep/LIS_output/LIS.CABLE.200001-200012.nc"
    wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_14Aug/WRF_output_copy/wrfout_d01_2013-01-01_03:00:00"
    var_names  = ["Rainf_tavg"] #var_LIS_soil_names # var_LIS_names # var_LIS_soil_names

    # for layer in np.arange(6):
    plot_map_var_lis(file_path, wrf_path, case_name, var_names, is_lnd=is_lnd)#, layer=layer)
