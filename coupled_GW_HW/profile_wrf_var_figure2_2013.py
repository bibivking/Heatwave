#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from netCDF4 import Dataset
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from scipy.interpolate import griddata
from wrf import (getvar, to_np, vertcross, CoordPair,
                 get_cartopy, latlon_coords, ALL_TIMES)
from common_utils import *

def get_time_cood(file_path, time_s ,time_e):

    # ****************** process time ******************
    ncfile      = Dataset(file_path)

    encoding    = 'utf-8' # Times in WRF output is btype, convert to string
    Time_s      = time_s - datetime(2000,1,1,0,0,0)
    Time_e      = time_e - datetime(2000,1,1,0,0,0)
    ntime       = len(ncfile.variables['Times'][:,0])
    time_tmp    = []

    # ****************** add three time coordiation ******************
    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(ncfile.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1,0,0,0))
    time = np.array(time_tmp)

    # for all days
    time_cood_all = (time>=Time_s) & (time<Time_e)
    time_t        = time[time_cood_all]
    doy_all       = [ time_t[i].days for i in np.arange(len(time_t)) ]
    time_t        = None
    
    # for daytime
    seconds       = [6.*60.*60.,18.*60.*60.]
    
    time_cood_day = []
    for j in np.arange(len(time)):
        if_day    = (time[j].seconds >= seconds[0]) & (time[j].seconds < seconds[1])
        time_cood_day.append((time[j]>=Time_s) & (time[j]<Time_e) & if_day)
    time_t        = time[time_cood_day]
    doy_day       = [ time_t[i].days for i in np.arange(len(time_t)) ]
    time_t        = None
    
    # for nighttime
    seconds       = [18.*60.*60.,6.*60.*60.]
    
    time_cood_night = []
    for j in np.arange(len(time)):
        if_night  = (time[j].seconds >= seconds[0]) | (time[j].seconds < seconds[1])
        time_cood_night.append((time[j]>=Time_s) & (time[j]<Time_e) & if_night)
    time_t        = time[time_cood_night]
    doy_night     = [ time_t[i].days for i in np.arange(len(time_t)) ]
    time_t        = None
    
    return time_cood_all, time_cood_day, time_cood_night, doy_all, doy_day, doy_night

def read_wrf_var(file_path):

    ncfile  = Dataset(file_path)
    Z       = getvar(ncfile, "z", timeidx=ALL_TIMES)
    Wa      = getvar(ncfile, "wa", units="m s-1", timeidx=ALL_TIMES)
    Ua      = getvar(ncfile, "ua", units="m s-1", timeidx=ALL_TIMES)
    T       = getvar(ncfile, 'temp', units='degC', timeidx=ALL_TIMES)
    S       = getvar(ncfile, 'QVAPOR', timeidx=ALL_TIMES)

    return Z, Wa, Ua, T, S

def get_time_masked(Z,T,S,Ua,Wa,time_cood):

    z  = Z[time_cood,:,:,:]
    t  = T[time_cood,:,:,:]
    s  = S[time_cood,:,:,:]
    ua = Ua[time_cood,:,:,:]
    wa = Wa[time_cood,:,:,:]

    return z, t, s, ua, wa

def get_average(Z,T,S,Ua,Wa,time_cood):

    z  = np.nanmean(Z[time_cood,:,:,:],axis=0)
    t  = np.nanmean(T[time_cood,:,:,:],axis=0)
    s  = np.nanmean(S[time_cood,:,:,:],axis=0)
    ua = np.nanmean(Ua[time_cood,:,:,:],axis=0)
    wa = np.nanmean(Wa[time_cood,:,:,:],axis=0)

    return z, t, s, ua, wa

def get_vertcross(file_path, z, t, s, ua, wa, lat_slt, lon_min, lon_max, doy, seconds=None):

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.

    ntime       = np.shape(z)[0]

    ncfile      = Dataset(file_path)
    start_point = CoordPair(lat=lat_slt, lon=lon_min)
    end_point   = CoordPair(lat=lat_slt, lon=lon_max)

    # ****************** calc interpolation ******************
    t_out     = np.zeros((ntime, 51, 48))
    s_out     = np.zeros((ntime, 51, 48))
    ua_out    = np.zeros((ntime, 51, 48))
    wa_out    = np.zeros((ntime, 51, 48))

    for i in np.arange(ntime):

        print("i = ", i)

        t_crs  = vertcross(t[i], z[i], wrfin=ncfile, start_point=start_point,
                            end_point=end_point, latlon=True, meta=True)
        s_crs  = vertcross(s[i], z[i], wrfin=ncfile, start_point=start_point,
                            end_point=end_point, latlon=True, meta=True)
        ua_crs = vertcross(ua[i], z[i], wrfin=ncfile, start_point=start_point,
                            end_point=end_point, latlon=True, meta=True)
        wa_crs = vertcross(wa[i], z[i], wrfin=ncfile, start_point=start_point,
                            end_point=end_point, latlon=True, meta=True)

        if i == 0:
            print("i == 0")
            loct = np.linspace(lon_min, lon_max, len(t_crs.coords['xy_loc']))
            print("loct", loct)
            vrt = np.arange(0, 5100., 100.)

        t_out[i], s_out[i], ua_out[i], wa_out[i] = \
                get_interpolation(t_crs, s_crs, ua_crs, wa_crs, loct, vrt)
        t_crs, s_crs, ua_crs, wa_crs = None, None, None, None
        print(loct)
        
    if seconds == None:
        t_cross  = np.nanmean(t_out, axis=0)
    elif seconds[0] < seconds[1]:
        # daytime - Tmax
        day_num = len(np.unique(doy))
        t       = np.zeros((day_num,len(t_out[0,:,0]),len(t_out[0,0,:])))
        for i in np.arange(day_num):
            is_the_day = [ doy[j] == np.unique(doy)[i] for j in np.arange(len(doy)) ]
            t[i,:,:]   = np.nanmax(t_out[is_the_day,:,:],axis=0)
        t_cross  = np.nanmean(t, axis=0)
    elif seconds[1] < seconds[0]:
        # night - Tmin
        day_num = len(np.unique(doy))
        t       = np.zeros((day_num,len(t_out[0,:,0]),len(t_out[0,0,:])))
        for i in np.arange(day_num):
            is_the_day = [ doy[j] == np.unique(doy)[i] for j in np.arange(len(doy)) ]
            t[i,:,:]   = np.nanmin(t_out[is_the_day,:,:],axis=0)
        t_cross  = np.nanmean(t, axis=0)
                
    s_cross  = np.nanmean(s_out, axis=0)
    ua_cross = np.nanmean(ua_out, axis=0)
    wa_cross = np.nanmean(wa_out, axis=0)
    print(loct)

    return t_cross, s_cross, ua_cross, wa_cross, loct, vrt

def get_interpolation(t_crs, s_crs, ua_crs, wa_crs, loct, vrt):

    print("get_interpolation")
    grid_X, grid_Y = np.meshgrid(loct,vrt)
    vertical_tmp   = to_np(t_crs.coords['vertical'])[:]

    grid_x, grid_y = np.meshgrid(loct,vertical_tmp)
    x              = np.reshape(grid_x,-1)
    y              = np.reshape(grid_y,-1)

    t_out  = griddata((x, y), np.reshape(to_np(t_crs),-1), (grid_X, grid_Y), method="linear")
    s_out  = griddata((x, y), np.reshape(to_np(s_crs),-1), (grid_X, grid_Y), method="linear")
    ua_out = griddata((x, y), np.reshape(to_np(ua_crs),-1), (grid_X, grid_Y), method="linear")
    wa_out = griddata((x, y), np.reshape(to_np(wa_crs),-1), (grid_X, grid_Y), method="linear")

    return t_out, s_out, ua_out, wa_out

def plot_profile_wrf_wind(file_paths, time_s, time_e, message=None, lat_slt=36, lon_min=130, lon_max=160):

    # ****************** Get time coordiation ******************
    time_cood_all, time_cood_day, time_cood_night, doy_all, doy_day, doy_night = \
        get_time_cood(file_paths[0], time_s, time_e)
    print("time_cood_day",time_cood_day)
    print("time_cood_day",time_cood_night)

    # ****************** Get the WRF variables ******************
    Z1,Wa1,Ua1,T1,S1 = read_wrf_var(file_paths[0])
    print("Z1", Z1)

    # ****************** Get time masked ******************
    z1_day, t1_day, s1_day, ua1_day, wa1_day = \
                       get_time_masked(Z1, T1, S1, Ua1, Wa1, time_cood_day)
    print("z1_day",z1_day)

    z1_night, t1_night, s1_night, ua1_night, wa1_night = \
                       get_time_masked(Z1, T1, S1, Ua1, Wa1, time_cood_night)

    # ****************** vertcross, interpolate and mean ******************
    seconds         = [6.*60.*60.,18.*60.*60.]
    t1_day_crs, s1_day_crs, ua1_day_crs, wa1_day_crs, xy_loc, vertical =\
        get_vertcross(file_paths[0], z1_day, t1_day, s1_day, ua1_day,
                      wa1_day, lat_slt, lon_min, lon_max, doy_day, seconds)

    seconds         = [18.*60.*60.,6.*60.*60.]
    t1_night_crs, s1_night_crs, ua1_night_crs, wa1_night_crs, xy_loc, vertical =\
        get_vertcross(file_paths[0], z1_night, t1_night, s1_night, ua1_night,
                      wa1_night, lat_slt, lon_min, lon_max, doy_night, seconds)

    if len(file_paths) > 1:
        # ****************** read second file ******************
        Z2,Wa2,Ua2,T2,S2 = read_wrf_var(file_paths[1])

        z2_day, t2_day, s2_day, ua2_day, wa2_day = \
                        get_time_masked(Z2,T2,S2,Ua2,Wa2,time_cood_day)

        z2_night, t2_night, s2_night, ua2_night, wa2_night = \
                        get_time_masked(Z2, T2, S2, Ua2, Wa2, time_cood_night)
                        
        seconds         = [6.*60.*60.,18.*60.*60.]
        t2_day_crs, s2_day_crs, ua2_day_crs, wa2_day_crs, xy_loc, vertical =\
            get_vertcross(file_paths[1], z2_day, t2_day, s2_day, ua2_day, wa2_day,
                          lat_slt, lon_min, lon_max, doy_day, seconds)

        seconds         = [18.*60.*60.,6.*60.*60.]
        t2_night_crs, s2_night_crs, ua2_night_crs, wa2_night_crs, xy_loc, vertical =\
            get_vertcross(file_paths[1], z2_night, t2_night, s2_night, ua2_night, 
                          wa2_night, lat_slt, lon_min, lon_max, doy_night, seconds)

        # ************** calc difference **************
        t_day_crs  = t2_day_crs  - t1_day_crs
        s_day_crs  = s2_day_crs  - s1_day_crs
        ua_day_crs = ua2_day_crs - ua1_day_crs
        wa_day_crs = wa2_day_crs - wa1_day_crs

        t_night_crs  = t2_night_crs  - t1_night_crs
        s_night_crs  = s2_night_crs  - s1_night_crs
        ua_night_crs = ua2_night_crs - ua1_night_crs
        wa_night_crs = wa2_night_crs - wa1_night_crs
        print("t_day_crs", t_day_crs)

    else:

        t_day_crs  = t1_day_crs
        s_day_crs  = s1_day_crs
        ua_day_crs = ua1_day_crs
        wa_day_crs = wa1_day_crs

        t_night_crs  = t1_night_crs
        s_night_crs  = s1_night_crs
        ua_night_crs = ua1_night_crs
        wa_night_crs = wa1_night_crs

    # ****************** plotting ******************
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[15,10],sharex=True, sharey=True, squeeze=True)
    plt.subplots_adjust(wspace=0.05, hspace=0.05)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color']     = almost_black
    plt.rcParams['xtick.color']     = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']      = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor']  = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')

    # Make the contour plot for var
    if len(file_paths) > 1:
        scale = 1.
    else:
        scale = 20.

    levels    = [-1.,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.]

    # Day temperature
    contour   = ax[0,0].contourf(xy_loc, vertical, t_day_crs, levels=levels, cmap=get_cmap("coolwarm"),extend='both')
    # cb_var    = fig.colorbar(contour, ax=ax[0,0])
    # cb_var.set_label('ΔT ($\mathregular{^o}$C)', loc='center') # rotation=270,
    # cb_var.ax[0,0].tick_params(labelsize=12)
    q         = ax[0,0].quiver(xy_loc[::3], vertical[::3], ua_day_crs[::3,::3],
                              wa_day_crs[::3,::3], angles='xy', scale_units='xy',
                              scale=scale, pivot='middle', color="white")
    # ax[0,0].quiverkey(q, X=0.90, Y=-0.05, U=scale, label=str(scale)+' m/s', labelpos='E')
    ax[0,0].text(0.02, 0.95, "(a) ΔT$_{max}$", transform=ax[0,0].transAxes, verticalalignment='top', bbox=props) # fontsize=14, 
    # ax[0,0].set_xlabel("Longitude", fontsize=12)
    ax[0,0].set_ylabel("Geopotential Height (m)")#, fontsize=12)

    # Night temperature
    contour   = ax[0,1].contourf(xy_loc, vertical, t_night_crs, levels=levels, cmap=get_cmap("coolwarm"),extend='both')
    cb_var    = fig.colorbar(contour, ax=ax[0], pad=0.01, orientation="vertical", aspect=20, shrink=0.88)
    cb_var.set_label('ΔT (${^o}$C)', loc='center') # rotation=270,
    # cb_var.ax[0,1].tick_params(labelsize=12)
    q         = ax[0,1].quiver(xy_loc[::3], vertical[::3], ua_night_crs[::3,::3],
                              wa_night_crs[::3,::3], angles='xy', scale_units='xy',
                              scale=scale, pivot='middle', color="white")
    # ax[0,1].quiverkey(q, X=0.90, Y=-0.05, U=scale, label=str(scale)+' m/s', labelpos='E')
    ax[0,1].text(0.02, 0.95, "(b) ΔT$_{min}$", transform=ax[0,1].transAxes, verticalalignment='top', bbox=props) # fontsize=14, 
    # ax[0,1].set_xlabel("Longitude", fontsize=12)
    # ax[0,1].set_ylabel("Geopotential Height (m)", fontsize=12)

    # Day specific humidity
    contour   = ax[1,0].contourf(xy_loc, vertical, s_day_crs*1000., levels=levels, cmap=get_cmap("coolwarm"),extend='both')
    # cb_var    = fig.colorbar(contour, ax=ax[1,0])
    # cb_var.set_label('ΔS (g kg$\mathregular{^-1}$)', loc='center') # rotation=270,
    # cb_var.ax[1,0].tick_params(labelsize=12)
    q         = ax[1,0].quiver(xy_loc[::3], vertical[::3], ua_day_crs[::3,::3],
                              wa_day_crs[::3,::3], angles='xy', scale_units='xy',
                              scale=scale, pivot='middle', color="white")
    # ax[1,0].quiverkey(q, X=0.90, Y=-0.05, U=scale, label=str(scale)+' m/s', labelpos='E')
    ax[1,0].text(0.02, 0.95, "(c) ΔS$_{day}$", transform=ax[1,0].transAxes, verticalalignment='top', bbox=props) # fontsize=14, 
    ax[1,0].set_xlabel("Longitude")#, fontsize=12)
    ax[1,0].set_ylabel("Geopotential Height (m)")#, fontsize=12)


    # Day specific humidity
    contour   = ax[1,1].contourf(xy_loc, vertical, s_night_crs*1000., levels=levels, cmap=get_cmap("coolwarm"),extend='both')
    cb_var    = fig.colorbar(contour, ax=ax[1], pad=0.01, orientation="vertical", aspect=20, shrink=0.88)
    cb_var.set_label('ΔS (g kg$^{-1}$)', loc='center') # rotation=270,
    # cb_var.ax[1,1].tick_params(labelsize=12)
    q         = ax[1,1].quiver(xy_loc[::3], vertical[::3], ua_night_crs[::3,::3],
                              wa_night_crs[::3,::3], angles='xy', scale_units='xy',
                              scale=scale, pivot='middle', color="white")
    ax[1,1].quiverkey(q, X=0.99, Y=-0.11, U=scale, label=str(scale)+' m/s', labelpos='E', color="black")
    ax[1,1].text(0.02, 0.95, "(d) ΔS$_{night}$", transform=ax[1,1].transAxes, verticalalignment='top', bbox=props) # fontsize=14, 
    ax[1,1].set_xlabel("Longitude")#, fontsize=12)
    # ax[1,1].set_ylabel("Geopotential Height (m)", fontsize=12)

    fig.savefig("./plots/figures/profile_wrf_Wind_"+message, bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    hw_name     = "hw2013_3Nov"
    lat_slt     = -36.
    lon_min     = 139.0
    lon_max     = 152.0

    # 2009
    if hw_name == "hw2009_3Nov":
        start_date= "20090122"
        end_date  = "20090213"
        time_s = datetime(2009,1,28,0,0,0,0)
        # time_e = datetime(2009,1,28,11,59,0,0)
        time_e = datetime(2009,2,8,23,59,0,0)
        # Time_s = datetime(2009,1,22,0,0,0,0)
        # Time_e = datetime(2009,2,13,23,59,0,0)

    elif hw_name == "hw2013_3Nov":
        start_date= "20121229"
        end_date  = "20130122"
        time_s = datetime(2013,1,4,0,0,0,0)
        time_e = datetime(2013,1,18,23,59,0,0)
        # Time_s = datetime(2012,12,29,0,0,0,0)
        # Time_e = datetime(2013,1,22,23,59,0,0)

    elif hw_name == "hw2019_3Nov":
        start_date= "20190108"
        end_date  = "20190130"
        time_s = datetime(2019,1,14,0,0,0)
        time_e = datetime(2019,1,26,23,59,0,0)
        # Time_s = datetime(2019,1,8,14,0,0,0)
        # Time_e = datetime(2019,1,30,0,0,0,0)
        #
    cpl_atmo_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+hw_name+'/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_'+start_date+'-'+end_date+'_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_'+start_date+'-'+end_date+'_fd'  # atmo output of wrf-cable run

    file_paths        = [cpl_atmo_file_fd,cpl_atmo_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw

    if len(file_paths) > 1:
        message = "GW-FD_"+str(time_s)+"-"+str(time_e)
    else:
        message = "GW_"+str(time_s)+"-"+str(time_e)

    plot_profile_wrf_wind(file_paths, time_s, time_e, message=message, lat_slt=lat_slt, lon_min=lon_min, lon_max=lon_max)
