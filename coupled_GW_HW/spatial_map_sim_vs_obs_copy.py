#!/usr/bin/python

from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_scale_offline
from common_utils import get_reverse_colormap

def tree_mask(file_path,pft_var_name):

    var = Dataset(file_path, mode='r')
    if pft_var_name == "Landcover_inst":
        pft = var.variables[pft_var_name][0,:,:]
    elif pft_var_name == "iveg":
        pft = var.variables[pft_var_name][:,:]

    return pft

def mask_by_lat_lon(file_path, loc_lat, loc_lon, lat_name, lon_name):

    '''
    make mask for the selected region
    '''

    file = nc.Dataset(file_path, mode='r')
    lat  = file.variables[lat_name][:]
    lon  = file.variables[lon_name][:]

    # print(lat)
    # print(lon)

    if len(np.shape(lat)) == 1:
        print("len(np.shape(lat)) == 1")
        lat_spc = lat[1] - lat[0]
        lon_spc = lon[1] - lon[0]
        lons, lats = np.meshgrid(lon, lat)
        mask  = (lats > (loc_lat[0] - lat_spc/2)) & (lats < (loc_lat[1] + lat_spc/2)) & (lons > (loc_lon[0] - lon_spc/2)) & (lons < (loc_lon[1] + lon_spc/2))
    elif len(np.shape(lat)) == 2:
        print("len(np.shape(lat)) == 2")
        ### caution: lat=100, lon=100 is a random pixel, lis run over a small domain may not have such a point
        lat_spc = lat[100,100] - lat[99,100]
        lon_spc = lon[100,100] - lon[100,99]
        print(lat_spc)
        print(lon_spc)
        ### caution: due to irregular space in lis, using lat/lon +lat/lon_spc/2 may includes more than 1 pixel.
        ### I therefore let the space divied by 2.1 rather than 2
        mask  = (lat > (loc_lat[0] - lat_spc/2.1)) & (lat < (loc_lat[1] + lat_spc/2.1)) & (lon > (loc_lon[0] - lon_spc/2.1)) & (lon < (loc_lon[1] + lon_spc/2.1))

    print(np.shape(mask))
    return mask

def read_obs_var(file_path, var_name, loc_lat=None, loc_lon=None, lat_name=None, lon_name=None):

    '''
    Read observation data, output time coordinate and variable array
    '''

    print(var_name)

    obs_file   = Dataset(file_path, mode='r')
    time_tmp   = nc.num2date(obs_file.variables['time'][:],obs_file.variables['time'].units,
                 only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time       = time_tmp - datetime(2000,1,1)
    ntime      = len(time)

    if loc_lat == None:
        Var_tmp = obs_file.variables[var_name][:]
        def_val = obs_file.variables[var_name]._FillValue
        Var = np.where(Var_tmp == def_val, np.nan, Var_tmp)
    else:
        # selected region
        if var_name == lat_name or var_name == lon_name:
            # read lat or lon
            mask = mask_by_lat_lon(file_path, loc_lat, loc_lon, lat_name, lon_name)
            lat  = obs_file.variables[lat_name]
            lon  = obs_file.variables[lon_name]
            if len(np.shape(lat)) == 1:
                lons, lats = np.meshgrid(lon, lat)
                if var_name == lat_name:
                    Var = np.where(mask,lats,np.nan)
                if var_name == lon_name:
                    Var = np.where(mask,lons,np.nan)
                print(np.shape(Var))
            elif len(np.shape(lat)) == 2:
                Var = np.where(mask, obs_file.variables[var_name][:], np.nan)
                print(np.shape(Var))
        else:
            # read var except lat or lat
            mask = mask_by_lat_lon(file_path, loc_lat, loc_lon, lat_name, lon_name)
            mask_multi = [ mask ] * ntime
            if var_name in ['E','Ei','Es','Et']:
                # change GLEAM's coordinates from (time, lon, lat) to (time, lat, lon)
                tmp = np.moveaxis(obs_file.variables[var_name], -1, 1)
            else:
                tmp = obs_file.variables[var_name][:]

            Var_tmp = np.where(mask_multi,tmp,np.nan)
            print(np.shape(Var_tmp))
            def_val = obs_file.variables[var_name]._FillValue
            Var = np.where(Var_tmp == def_val, np.nan, Var_tmp)

    return time,Var

def read_off_var(file_path, var_name, loc_lat=None, loc_lon=None, lat_name=None, lon_name=None):

    '''
    Read observation data, output time coordinate and variable array
    '''

    print(var_name)

    off_file   = Dataset(file_path, mode='r')
    time_tmp   = nc.num2date(off_file.variables['time'][:],off_file.variables['time'].units,
                 only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time       = time_tmp - datetime(2000,1,1)
    ntime      = len(time)

    if loc_lat == None:
        Var_tmp = off_file.variables[var_name][:]
        def_val = off_file.variables[var_name]._FillValue
        Var = np.where(Var_tmp == def_val, np.nan, Var_tmp)
    else:
        # selected region
        if var_name == lat_name or var_name == lon_name:
            # read lat or lon
            mask = mask_by_lat_lon(file_path, loc_lat, loc_lon, lat_name, lon_name)
            lat  = off_file.variables[lat_name]
            lon  = off_file.variables[lon_name]
            if len(np.shape(lat)) == 1:
                lons, lats = np.meshgrid(lon, lat)
                if var_name == lat_name:
                    Var = np.where(mask,lats,np.nan)
                if var_name == lon_name:
                    Var = np.where(mask,lons,np.nan)
                print(np.shape(Var))
            elif len(np.shape(lat)) == 2:
                Var_tmp = np.where(mask, off_file.variables[var_name][:], np.nan)
                def_val = off_file.variables[var_name]._FillValue
                Var = np.where(Var_tmp == def_val, np.nan, Var_tmp)
                print(np.shape(Var))
        else:
            # read var except lat or lat
            mask = mask_by_lat_lon(file_path, loc_lat, loc_lon, lat_name, lon_name)
            mask_multi = [ mask ] * ntime
            if var_name in ['E','Ei','Es','Et']:
                # change GLEAM's coordinates from (time, lon, lat) to (time, lat, lon)
                tmp = np.moveaxis(off_file.variables[var_name], -1, 1)
            else:
                tmp = off_file.variables[var_name][:]

            Var_tmp = np.where(mask_multi,tmp,np.nan)
            print(np.shape(Var_tmp))
            def_val = off_file.variables[var_name]._FillValue
            Var = np.where(Var_tmp == def_val, np.nan, Var_tmp)

    return time,Var

def spital_var(time,Var,year_s,year_e):

    Year_s = datetime(year_s,1,1)   - datetime(2000,1,1)
    Year_e = datetime(year_e,12,31) - datetime(2000,1,1)
    time_cood = (time>=Year_s) & (time<=Year_e)
    print(np.shape(Var))
    var = np.nanmean(Var[time_cood,:,:],axis=0)

    return var

def time_series_var(time,Var,year_s,year_e):

    Year_s = datetime(year_s,1,1)   - datetime(2000,1,1)
    Year_e = datetime(year_e,12,31) - datetime(2000,1,1)

    var_tmp  = Var[(time>=Year_s) & (time<=Year_e),:,:]
    var      = np.nanmean(var_tmp,axis=(1,2))
    Time     = time[(time>=Year_s) & (time<=Year_e)]

    return Time,var

def plot_spital_map(file_path, var_name, year_s, year_e, loc_lat=None, loc_lon=None, lat_name=None, lon_name=None, message=None):

    print("======== In plot_spital_map =========")

    # Open the NetCDF4 file (add a directory path if necessary) for reading:

    time, Var  = read_obs_var(file_path, var_name, loc_lat, loc_lon, lat_name, lon_name)
    time, lats = read_obs_var(file_path, lat_name, loc_lat, loc_lon, lat_name, lon_name)
    time, lons = read_obs_var(file_path, lon_name, loc_lat, loc_lon, lat_name, lon_name)

    var       = spital_var(time,Var,year_s,year_e)
    print(var)

    fig = plt.figure(figsize=(7,5))
    ax = plt.axes(projection=ccrs.PlateCarree())

    if loc_lat == None:
        ax.set_extent([140,154,-40,-28])
    else:
        ax.set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])

    ax.coastlines(resolution="50m",linewidth=1)
    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = True

    if loc_lat == None:
        gl.xlocator = mticker.FixedLocator([140,145,150])
        gl.ylocator = mticker.FixedLocator([-40,-35,-30])
    else:
        gl.xlocator = mticker.FixedLocator(loc_lon)
        gl.ylocator = mticker.FixedLocator(loc_lat)

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':'black'}
    gl.ylabel_style = {'size':10, 'color':'black'}
    # Plot windspeed

    if var_name in ['Qs','Qsb','Rainf','Evap','ESoil','ECanop','TVeg']:
        scale = 24.*3600.*365
    else:
        scale = 1.

    # clevs = np.linspace( np.min(var),np.max(var), num=20)
    plt.contourf(lons, lats, var*scale,  transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG) # clevs,
    plt.title(var_name, size=16)
    cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    # cb.set_label(units,size=14,rotation=270,labelpad=15)
    cb.ax.tick_params(labelsize=10)
    if message == None:
        message = var_name
    else:
        message = message + "_" + var_name
    plt.savefig('./plots/Oct2021/spatial_map_obs_'+message+'.png',dpi=300)

def plot_spital_map_diff(file_path1, file_path2, var_name, year_s, year_e,
                         loc_lat=None, loc_lon=None, lat_name=None, lon_name=None,
                         message=None):

    print("======== In plot_spital_map =========")

    # Open the NetCDF4 file (add a directory path if necessary) for reading:

    time1, Var1 = read_obs_var(file_path1, var_name, loc_lat, loc_lon, lat_name, lon_name)
    time2, Var2 = read_obs_var(file_path2, var_name, loc_lat, loc_lon, lat_name, lon_name)
    time1, lats1 = read_obs_var(file_path1, lat_name, loc_lat, loc_lon, lat_name, lon_name)
    time1, lons1 = read_obs_var(file_path1, lon_name, loc_lat, loc_lon, lat_name, lon_name)

    var1       = spital_var(time1,Var1,year_s,year_e)
    var2       = spital_var(time1,Var2,year_s,year_e)
    var        = var2-var1
    print(var)

    fig = plt.figure(figsize=(7,5))
    ax = plt.axes(projection=ccrs.PlateCarree())

    if loc_lat == None:
        ax.set_extent([140,154,-40,-28])
    else:
        ax.set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])

    ax.coastlines(resolution="50m",linewidth=1)
    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = True

    if loc_lat == None:
        gl.xlocator = mticker.FixedLocator([140,145,150])
        gl.ylocator = mticker.FixedLocator([-40,-35,-30])
    else:
        gl.xlocator = mticker.FixedLocator(loc_lon)
        gl.ylocator = mticker.FixedLocator(loc_lat)

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':'black'}
    gl.ylabel_style = {'size':10, 'color':'black'}
    # Plot windspeed

    if var_name in ['Qs','Qsb','Rainf','Evap','ESoil','ECanop','TVeg']:
        scale = 24.*3600.*365
    else:
        scale = 1.

    # clevs = np.linspace( np.min(var),np.max(var), num=20)
    plt.contourf(lons, lats, var*scale,  transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG) # clevs,
    plt.title(var_name, size=16)
    cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    # cb.set_label(units,size=14,rotation=270,labelpad=15)
    cb.ax.tick_params(labelsize=10)
    if message == None:
        message = var_name
    else:
        message = message + "_" + var_name
    plt.savefig('./plots/Oct2021/spatial_map_obs_'+message+'.png',dpi=300)

def plot_time_series(file_paths, var_names, year_s, year_e, loc_lat=None, loc_lon=None, lat_name=None, lon_name=None):

    print("======== In plot_time_series =========")

    fig, ax = plt.subplots()
    Time1, Var1 = read_obs_var(file_paths[0], var_names[0], loc_lat, loc_lon, lat_name[0], lon_name[0])
    time1, var1 = time_series_var(Time1,Var1,year_s,year_e)
    t1 = []
    for i in np.arange(len(time1)):
        t1.append(time1[i].days)
    if var_names[0] in ['Qs','Qsb','Rainf','Evap','ESoil','ECanop','TVeg']:
        scale = 24.*3600.
    else:
        scale = 1.

    # ax.plot(t1, var1*scale, c = "blue", label="var1", alpha=0.5)
    # if len(file_paths) > 1:
        # Time2, Var2 = read_obs_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_name[1], lon_name[1])
        # time2, var2 = time_series_var(Time2,Var2,year_s,year_e)
        # t2 = []
        # for i in np.arange(len(time2)):
        #     t2.append(time2[i].days)
        # if var_names[1] in ['Qs','Qsb','Rainf','Evap','ESoil','ECanop','TVeg']:
        #     scale = 24.*3600.
        # else:
        #     scale = 1.
        # ax.plot(t2, var2*scale, c = "red", label="var2", alpha=0.5)

    Time2, Var2 = read_obs_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_name[1], lon_name[1])
    time2, var2 = time_series_var(Time2,Var2,year_s,year_e)
    var = [[var1],[var2]]
    ax.plot(t1, var*scale, alpha=0.5)
    # ax.set_ylabel('mm')
    ax.set_title(var_names[0])
    # ax.set_xticks(x1[::1440])
    # ax.set_xticklabels(np.arange(year_s,year_e,1))
    ax.legend()

    fig.tight_layout()
    message = var_names[0]
    if loc_lat != None:
        message = message + "_lat="+str(loc_lat[0]) +"-"+str(loc_lat[1]) + "_lon="+str(loc_lon[0])+"-"+str(loc_lon[1])

    plt.savefig('./plots/Oct2021/time_series_lis_vs_off_'+message+'.png',dpi=300)

if __name__ == "__main__":

    # #######################
    #   plot_spital_map     #
    # #######################

    DOLCE_path = "/g/data/w35/mm3972/data/DOLCE/v3/"
    DOLCE_file = DOLCE_path+"DOLCE_v3_2000-2018.nc"
    GLEAM_path = "/g/data/w35/Shared_data/Observations/Global_ET_products/GLEAM_v3_5/v3-5a/yearly/"
    GLEAM_file = GLEAM_path + "E_1980-2020_GLEAM_v3.5a_YR.nc"
    GRACE_path = "/g/data/w35/mm3972/data/GRACE/GRACE_JPL_RL06/GRACE_JPLRL06M_MASCON/"
    GRACE_file = GRACE_path + "GRCTellus.JPL.200204_202004.GLO.RL06M.MSCNv02CRI.nc"

    GW_off_path   = "/g/data/w35/mm3972/model/cable/runs/sanaa_run/test_fix_satfrac_10km/outputs/"
    GW_off_file   = GW_off_path + "cable_out_2000-2019.nc"

    FD_off_path   = "/g/data/w35/mm3972/model/cable/runs/sanaa_run/test_fix_satfrac_10km_fd/outputs/"
    FD_off_file   = FD_off_path + "cable_out_2000-2019.nc"

    year_s     = 2000
    year_e     = 2019
    loc_lat    = [-40,-28]
    loc_lon    = [140,154]

    var_names   = ['Fwsoil','Qs','Qsb','WatTable','Qle','Qh','Qg','GWMoist','Rainf','Evap','ESoil','ECanop','TVeg']
    lat_name   = "latitude"#"lat"
    lon_name   = "longitude"#"lon"
    for var_name in var_names:
        plot_spital_map_diff(FD_off_file, GW_off_file, var_name, year_s, year_e,
                         loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name,
                         message="GW-FD")
        # plot_spital_map( FD_off_file, var_name, year_s, year_e,
        #                  loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name,
        #                  message="FD")

    # var_name   = "E" #"E"#'lwe_thickness'#'E'
    # lat_name   = "lat"
    # lon_name   = "lon"
    # plot_spital_map(GLEAM_file, var_name, year_s, year_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name)
    #
    # file_path  = [GW_off_file, FD_off_file]
    # lat_name   = ["latitude","latitude"]
    # lon_name   = ["longitude","longitude"]
    #
    # var_names  = [ ["GWMoist","GWMoist"],
    #                ["Evap","Evap"],
    #                ["TVeg","TVeg"],
    #                ["ESoil","ESoil"],
    #                ["ECanop","ECanop"],
    #                ["Qs","Qs"],
    #                ["Qsb","Qsb"],
    #                ["Rainf","Rainf"],
    #                ["WatTable","WatTable"],
    #                ["Qle","Qle"],
    #                ["Qh","Qh"],
    #                ["Qg","Qg"],
    #                ["RadT","RadT"],
    #                ["VegT","VegT"],
    #                ["Fwsoil","Fwsoil"]]
    #
    # for var_num in np.arange(len(var_names)):
    #     plot_time_series(file_path, var_names[var_num], year_s, year_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name)
