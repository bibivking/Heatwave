#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale

def plot_map_landinfo_ts(file_paths, wrf_path, case_names, var_names, ts, layer=None):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var1 = Dataset(file_paths[0], mode='r')
    wrf = Dataset(wrf_path,  mode='r')

    # use WRF output's lat & lon, since LIS output has default value
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    for var_name in var_names:
        print(var_name)
        scale, units = get_land_var_scale(var_name)
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

        cmap = plt.cm.jet
        if var_name in ["Landcover_inst","Soiltype_inst"]:
            clevs = np.linspace(0.5,15.5, num=16)
            plt.contourf(lon, lat, Var[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        elif var_name in ["Albedo_inst","SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst"]:
            clevs = np.linspace(0.0,0.5, num=11)
            plt.contourf(lon, lat, Var[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        elif var_name in ["SandFrac_inst","ClayFrac_inst","SiltFrac_inst"]:
            clevs = np.linspace(0.0,1., num=11)
            plt.contourf(lon, lat, Var[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        elif var_name in ["LAI_inst"]:
            clevs = np.linspace(0.,8., num=16)
            plt.contourf(lon, lat, Var[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        else:
            plt.contourf(lon, lat, Var[:,:], transform=ccrs.PlateCarree(),cmap=cmap)

        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

        if units == None:
            units_string = var1.variables[var_name].units
        else:
            units_string = units

        plt.title(var_name+' ('+units_string+') ts= '+str(ts), size=16)

        cb.ax.tick_params(labelsize=10)
        cb.set_label(units_string,size=14,rotation=270,labelpad=15)

        if len(np.shape(Var)) == 3:
            plt.savefig('./plots/spatial_map_landinfo_'+var_name+'_lvl-'+str(layer)+'_'+case_names[0]+'.png',dpi=300)
        else:
            plt.savefig('./plots/spatial_map_landinfo_'+var_name+'_'+case_names[0]+'.png',dpi=300)
        Var = None


def plot_map_landinfo_offline_CABLE(file_path, var_names):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var = Dataset(file_path, mode='r')

    lons = var.variables['longitude']
    lats = var.variables['latitude']
    lon, lat = np.meshgrid(lons, lats)

    for var_name in var_names:
        print(var_name)
        Var = var.variables[var_name]

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

        cmap = plt.cm.jet
        if var_name in ["iveg","isoil"]:
            clevs = np.linspace(0.5,15.5, num=16)
            plt.contourf(lon, lat, Var[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        elif var_name in ["Albedo","sfc","ssat","swilt"]:
            clevs = np.linspace(0.0,0.5, num=11)
            plt.contourf(lon, lat, Var[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        elif var_name in ["sand","clay","silt"]:
            clevs = np.linspace(0.0,1., num=11)
            plt.contourf(lon, lat, Var[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        elif var_name in ["LAI"]:
            clevs = np.linspace(0.,8., num=16)
            plt.contourf(lon, lat, Var[11,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        elif var_name in ["sucs","hyds"]:
            plt.contourf(lon, lat, Var[:,:]*1000., transform=ccrs.PlateCarree(),cmap=cmap)
        else:
            plt.contourf(lon, lat, Var[:,:], transform=ccrs.PlateCarree(),cmap=cmap)

        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

        # units_string = var.variables[var_name].units


        plt.title(var_name, size=16)

        cb.ax.tick_params(labelsize=10)
        # cb.set_label(units_string,size=14,rotation=270,labelpad=15)


        plt.savefig('./plots/spatial_map_landinfo_offline_cable_'+var_name+'_zobler.png',dpi=300)

        Var = None

if __name__ == "__main__":

    ##################################
    #   plot plot_map_landinfo_ts    # 
    ##################################
    var_landinfo_3D_names =  ["Landcover_inst","Soiltype_inst","SandFrac_inst","ClayFrac_inst","SiltFrac_inst",
                              "SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst","Hyds_inst","Bch_inst","Sucs_inst",
                              "Elevation_inst","LAI_inst","Albedo_inst"]

    case_names = ['ctl_14Aug']
    file_name  = "LIS.CABLE.198212-201301.nc"
    file_paths = []
    is_diff    = False #False
    for case_name in case_names:
        path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
        file_path  = path + file_name
        file_paths.append(file_path)
    
    case_name = 'ctl_25Jul'
    # Since lon and lat in LIS contain default values, to use plt.contourf, I take lon/lat from WRF output
    wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/wrfout_d01_2013-01-01_03:00:00"

    d2012121   = 10524
    ts         = d2012121

    var_names  = var_landinfo_3D_names
    plot_map_landinfo_ts(file_paths, wrf_path, case_names, var_names, ts)

    #############################################
    #   plot plot_map_landinfo_offline_CABLE    # 
    #############################################
    file_path = "/g/data/w35/mm3972/model/cable/src/CABLE-AUX/offline/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc"
    var_names  = ["iveg","isoil","sand","clay","silt",
                  "sfc","ssat","swilt","hyds","bch","sucs",
                  "LAI","Albedo",'elevation']
    plot_map_landinfo_offline_CABLE(file_path,var_names)
