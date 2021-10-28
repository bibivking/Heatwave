#!/usr/bin/python

from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
from wrf import (getvar, interplevel,ALL_TIMES)

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
        # print("loc 1")
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
        # print("loc 2")
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
            # print("loc 3")
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

def read_wrf_var(file_path, var_name, var_unit, height, loc_lat=None, loc_lon=None):

    print("read "+var_name+" from wrf output")
    wrf_file = Dataset(file_path)

    Var      = getvar(wrf_file, var_name, units=var_unit,timeidx=ALL_TIMES)
    p        = getvar(wrf_file, "pressure",timeidx=ALL_TIMES)
    var_tmp  = interplevel(Var, p, height)
    ntime    = len(p[:,0,0,0])
    if loc_lat == None:
        var  = var_tmp
    else:
        mask = mask_by_lat_lon(file_path, loc_lat, loc_lon, 'XLAT', 'XLONG')
        mask_multi = [ mask ] * ntime
        var = np.where(mask_multi,var_tmp,np.nan)

    return var

# ========================= Spitial & temporal Average =========================
def spital_var(time,Var,time_s,time_e):

    Time_s = time_s - datetime(2000,1,1,0,0,0)
    Time_e = time_e - datetime(2000,1,1,0,0,0)
    time_cood = (time>=Time_s) & (time<Time_e)

    # print('===== np.shape(Var[time_cood,:,:]) =====')
    # print(np.shape(Var[time_cood,:,:]))

    var = np.nanmean(Var[time_cood,:,:],axis=0)
    # for rain
    # var = np.sum(Var[time_cood,:,:],axis=0)
    # print('===== np.shape(var) =====')
    # print(np.shape(var))

    # np.savetxt("test_var.txt",var,delimiter=",")
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
                        "SnowDepth_inst","SoilWet_inst",
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
