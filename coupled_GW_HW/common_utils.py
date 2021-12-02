#!/usr/bin/python

from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
import matplotlib.colors as colors
from datetime import datetime, timedelta
from wrf import (getvar, interplevel, get_cartopy, cartopy_xlim,
                 cartopy_ylim, to_np, latlon_coords, ALL_TIMES)

# =============================== Basci Functions ==============================
def leap_year(year):
   if (year % 4) == 0:
       return 366
   else:
       return 365

def get_scale(var_name):
    pa2hpa     = ['ps']
    m2mm       = ['tp']
    mm_s2mm_yr = ['Qs','Qsb','Rainf','Evap','ESoil','ECanop','TVeg']
    mm_s2mm_day= ['Rainf_tavg']
    if var_name in mm_s2mm_yr:
        scale = 24.*3600.*365
    elif var_name in mm_s2mm_day:
        scale = 24.*3600.
    elif var_name in pa2hpa:
        scale = 0.01
    elif var_name in m2mm:
        scale = 1000.
    else:
        scale = 1.
    return scale

def UTC_to_AEST(time):

    Time = time + timedelta(hours=10)

    return Time

# ===================================== Mask ===================================
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
    if len(np.shape(file.variables[lat_name][:])) == 3:
        # print("len(np.shape(file.variables[lat_name][:])) == 3")
        # print(lat_name)
        lat  = file.variables[lat_name][0,:,:]
        lon  = file.variables[lon_name][0,:,:]
    else:
        lat  = file.variables[lat_name][:]
        lon  = file.variables[lon_name][:]

    # print(lat)
    # print(lon)

    if len(np.shape(lat)) == 1:
        # print("len(np.shape(lat)) == 1")
        lat_spc = lat[1] - lat[0]
        lon_spc = lon[1] - lon[0]
        lons, lats = np.meshgrid(lon, lat)
        mask  = (lats > (loc_lat[0] - lat_spc/2)) & (lats < (loc_lat[1] + lat_spc/2)) & (lons > (loc_lon[0] - lon_spc/2)) & (lons < (loc_lon[1] + lon_spc/2))
    elif len(np.shape(lat)) == 2:
        # print("len(np.shape(lat)) == 2")
        ### caution: lat=100, lon=100 is a random pixel, lis run over a small domain may not have such a point
        lat_spc = lat[50,50] - lat[49,50]
        lon_spc = lon[50,50] - lon[50,49]
        # print(lat_spc)
        # print(lon_spc)
        ### caution: due to irregular space in lis, using lat/lon +lat/lon_spc/2 may includes more than 1 pixel.
        ### I therefore let the space divied by 2.1 rather than 2
        mask  = (lat > (loc_lat[0] - lat_spc/2.1)) & (lat < (loc_lat[1] + lat_spc/2.1)) & (lon > (loc_lon[0] - lon_spc/2.1)) & (lon < (loc_lon[1] + lon_spc/2.1))

    # print(np.shape(mask))
    return mask

def time_mask(time, time_s, time_e, seconds=None):

    Time_s = time_s - datetime(2000,1,1,0,0,0)
    Time_e = time_e - datetime(2000,1,1,0,0,0)

    if seconds == None:
        time_cood = (time>=Time_s) & (time<Time_e)
    else:
        time_cood = []
        for j in np.arange(len(time)):
            if seconds[0] >= seconds[1]:
                if_seconds = (time[j].seconds >= seconds[0]) | (time[j].seconds < seconds[1])
            else:
                if_seconds = (time[j].seconds >= seconds[0]) & (time[j].seconds < seconds[1])
            time_cood.append( (time[j]>=Time_s) & (time[j]<Time_e) & if_seconds)

    return time_cood

# ================================ Read variables ==============================
def read_var(file_path, var_name, loc_lat=None, loc_lon=None, lat_name=None, lon_name=None):

    '''
    Read observation data, output time coordinate and variable array
    '''

    print(var_name)

    obs_file   = Dataset(file_path, mode='r')
    time_tmp   = nc.num2date(obs_file.variables['time'][:],obs_file.variables['time'].units,
                 only_use_cftime_datetimes=False, only_use_python_datetimes=True)

    time       = time_tmp - datetime(2000,1,1,0,0,0)

    ntime      = len(time)

    if loc_lat == None:
        Var_tmp = obs_file.variables[var_name][:]
        if hasattr(obs_file.variables[var_name], '_FillValue'):
            # hasattr(a,"b"): check whether object a has attribute 'b'
            def_val = obs_file.variables[var_name]._FillValue
            Var = np.where(Var_tmp == def_val, np.nan, Var_tmp)
        elif hasattr(obs_file.variables[var_name], '_fillvalue'):
            def_val = obs_file.variables[var_name]._fillvalue
            Var = np.where(Var_tmp == def_val, np.nan, Var_tmp)
        else:
            Var = Var_tmp
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
                # print(np.shape(Var))
            elif len(np.shape(lat)) == 2:
                Var = np.where(mask, obs_file.variables[var_name][:], np.nan)
                # print(np.shape(Var))
            elif len(np.shape(lat)) == 3:
                Var = np.where(mask, obs_file.variables[var_name][0,:,:], np.nan)
                # print(np.shape(Var))
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
            # print(np.shape(Var_tmp))
            if hasattr(obs_file.variables[var_name], '_FillValue'):
                def_val = obs_file.variables[var_name]._FillValue
                Var = np.where(Var_tmp == def_val, np.nan, Var_tmp)
            elif hasattr(obs_file.variables[var_name], '_fillvalue'):
                def_val = obs_file.variables[var_name]._fillvalue
                Var = np.where(Var_tmp == def_val, np.nan, Var_tmp)
            else:
                Var = Var_tmp
    return time,Var

def read_wrf_surf_var(file_path, var_name, loc_lat=None, loc_lon=None, mask_map=None):

    # output: [time,lat,lon]
    print("read "+var_name+" from wrf output")

    var_3D = [
                'rh2',  # 2m Relative Humidity
                'T2',   # 2m Temperature
                'td2',  # 2m Dew Point Temperature
                'slp',  # Sea Level Pressure
                'ter',  # Model Terrain Height
                'ctt',  # Cloud Top Temperature
                'mdbz', # Maximum Reflectivity
                'pw',   # Precipitable Water
                'updraft_helicity', # Updraft Helicity
                'helicity',        # Storm Relative Helicity
              ]

    wrf_file = Dataset(file_path)
    p        = getvar(wrf_file, "pressure",timeidx=ALL_TIMES)
    if var_name == 'cape_2d':
        # 'cape_2d', # 2D CAPE (MCAPE/MCIN/LCL/LFC)
        var_tmp  = getvar(wrf_file, var_name, timeidx=ALL_TIMES)[0]
        print("======= cape_2d =======")
        print(var_tmp)
    elif var_name == 'cloudfrac':
        # 'cloudfrac', # Cloud Fraction
        var_tmp  = getvar(wrf_file, var_name, timeidx=ALL_TIMES)[0]
        print("======= cloudfrac =======")
        print(var_tmp)
    else:
        var_tmp  = getvar(wrf_file, var_name, timeidx=ALL_TIMES)

    if loc_lat == None:
        if mask_map is not None:
            ntime  = len(p[:,0,0,0])
            mask_multi = [ mask_map ] * ntime
            var  = np.where(mask_multi,var_tmp,np.nan)
        else:
            var  = var_tmp
    else:
        ### need to fix, not work
        ntime  = len(p[:,0,0,0])
        mask   = mask_by_lat_lon(file_path, loc_lat, loc_lon, 'XLAT', 'XLONG')
        if mask_map is not None:
            mask = np.where(np.all([mask,mask_map],axis=0), True, False)
        # np.savetxt("check",mask)
        mask_multi = [ mask ] * ntime
        var    = np.where(mask_multi,var_tmp,np.nan)

    return var

def read_wrf_hgt_var(file_path, var_name, var_unit=None, height=None, loc_lat=None, loc_lon=None):

    print("read "+var_name+" from wrf output")

    var_4D =  [
                'p',    # Full Model Pressure
                'avo',    # Absolute Vorticity
                'eth',    # Equivalent Potential Temperature
                'dbz',    # Reflectivity
                'geopt',  # Geopotential for the Mass Grid  
                'omg',  # Omega
                'pvo',  # Potential Vorticity
                'rh',   # Relative Humidity
                'td',   # Dew Point Temperature
                'tc',   # Temperature in Celsius
                'th',   # Potential Temperature
                'temp', # Temperature (in specified units)
                'tv',   # Virtual Temperature
                'twb',  # Wet Bulb Temperature
                'ua',   # U-component of Wind on Mass Points
                'va',   # V-component of Wind on Mass Points
                'wa',   # W-component of Wind on Mass Points
                'z',    # Model Height for Mass Grid
                'cape_3d',# 3D CAPE and CIN
                'height_agl', # Model Height for Mass Grid (AGL)
                ]


    wrf_file = Dataset(file_path)
    p        = getvar(wrf_file, "pressure",timeidx=ALL_TIMES)

    # if var_name in var_4D:
    if var_unit == None:
        if var_name == 'cape_3d':
            Var  = getvar(wrf_file, var_name, timeidx=ALL_TIMES)[0]
            print("======= Var =======")
            print(Var)
        else:
            Var  = getvar(wrf_file, var_name, timeidx=ALL_TIMES)
    else:
        Var      = getvar(wrf_file, var_name, units=var_unit, timeidx=ALL_TIMES)
    # else:
    #     Var  = wrf_file.variables[var_name][:]

    if height == None:
        var_tmp  = Var
    else:
        var_tmp  = interplevel(Var, p, height) ## need to check whether works with wrf_file.variables[var_name][:]

    if loc_lat == None:
        var  = var_tmp
    else:
        # here only suit 2D and 3D var
        ### need to fix, not work
        ntime  = len(p[:,0,0,0])
        mask   = mask_by_lat_lon(file_path, loc_lat, loc_lon, 'XLAT', 'XLONG')
        # np.savetxt("check",mask)
        mask_multi = [ mask ] * ntime
        var    = np.where(mask_multi,var_tmp,np.nan)

    return var

# ========================= Spitial & temporal Average =========================
def spital_var(time, Var, time_s, time_e, seconds=None):

    time_cood = time_mask(time, time_s, time_e, seconds)
    var       = np.nanmean(Var[time_cood],axis=0)

    # np.savetxt("test_var.txt",var,delimiter=",")
    return var

def spital_var_max(time, Var, time_s, time_e, seconds=None):

    time_cood = time_mask(time, time_s, time_e, seconds)
    time_slt  = time[time_cood]

    var_slt  = Var[time_cood,:,:]

    days     = []
    var_tmp  = []

    for t in time_slt:
        days.append(t.days)

    for d in np.arange(days[0],days[-1]+1):
        var_tmp.append(np.nanmax(var_slt[days == d,:,:],axis=0))

    var = np.nanmean(var_tmp,axis=0)

    return var

def spital_var_min(time, Var, time_s, time_e, seconds=None):
    
    time_cood = time_mask(time, time_s, time_e, seconds)
    time_slt  = time[time_cood]

    var_slt  = Var[time_cood,:,:]

    days     = []
    var_tmp  = []

    for t in time_slt:
        days.append(t.days)

    for d in np.arange(days[0],days[-1]+1):
        var_tmp.append(np.nanmin(var_slt[days == d,:,:],axis=0))

    var = np.nanmean(var_tmp,axis=0)

    return var

def spital_ERAI_tp(time,Var,time_s,time_e):

    Time_s = time_s - datetime(2000,1,1,0,0,0)
    Time_e = time_e - datetime(2000,1,1,0,0,0)

    is_12_24 = []
    for i in np.arange(len(time)):
        is_12_24.append((time[i].seconds == 43200) | (time[i].seconds == 0))
    time_cood = (time>=Time_s) & (time<Time_e) & is_12_24
    # print(float(np.count_nonzero(time_cood==True)))
    var = np.nanmean(Var[time_cood,:,:],axis=0) * float(np.count_nonzero(time_cood==True))

    # print(np.nansum(var))
    return var

def time_series_var(time,Var,time_s,time_e):

    Time_s = time_s - datetime(2000,1,1,0,0,0)
    Time_e = time_e - datetime(2000,1,1,0,0,0)

    var_tmp  = Var[(time>=Time_s) & (time<=Time_e),:,:]
    var      = np.nanmean(var_tmp,axis=(1,2))
    Time     = time[(time>=Time_s) & (time<=Time_e)]

    return Time,var

def time_series_statistic(time,Var,time_s,time_e):

    Time_s = time_s - datetime(2000,1,1,0,0,0)
    Time_e = time_e - datetime(2000,1,1,0,0,0)

    var_tmp  = Var[(time>=Time_s) & (time<=Time_e),:,:]

    ##### how to keep these multi lines rather than average????
    var_mean = np.nanmean(var_tmp,axis=(1,2))
    var_max  = np.nanmax(var_tmp,axis=(1,2))
    var_min  = np.nanmin(var_tmp,axis=(1,2))

    return var_mean, var_max, var_min

def time_series_time(time,time_s,time_e):

    Time_s = time_s - datetime(2000,1,1,0,0,0)
    Time_e = time_e - datetime(2000,1,1,0,0,0)
    Time     = time[(time>=Time_s) & (time<=Time_e)]

    return Time

# =============================== Plots setting ================================

def get_reverse_colormap(var_name):

    '''
    To tell whether it needs a reversed colormap
    '''

    var_reverse_yes = [ "Rainf_f_inst","Rainf_tavg","Evap_tavg","ECanop_tavg","TVeg_tavg","ESoil_tavg",
                        "Qs_tavg","Qsb_tavg", "Snowf_tavg","GPP_tavg","Qle_tavg","SoilMoist_inst",
                        "FWsoil_tavg","SnowCover_inst","Qair_f_inst","Wind_f_inst","SWE_inst",
                        "SnowDepth_inst","SoilWet_inst", "EF",
                        "rh2"]
    var_reverse_no  = [ "Qh_tavg","Qg_tavg","Swnet_tavg","Lwnet_tavg","SWdown_f_inst","LWdown_f_inst",
                        "VegT_tavg","AvgSurfT_tavg","Tair_f_inst","SoilTemp_inst","Albedo_inst",
                        "Psurf_f_inst",
                        "T2"]

    if var_name in var_reverse_yes:
        return_value = True
    elif var_name in var_reverse_no:
        return_value = False
    else:
        return_value = None

    return return_value

def get_wrf_var_range_diff(var_name):

    '''
    Get range
    '''

    var_degc       = ["T2"]
    var_percent    = ["rh2"]
    ranges         = [0.0,0.0]

    if var_name in var_degc:
        ranges[0] = -2.
        ranges[1] = 2.
    elif var_name in var_percent:
        ranges[0] = -20.
        ranges[1] = 20.
    else:
        ranges = None
    return ranges

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
