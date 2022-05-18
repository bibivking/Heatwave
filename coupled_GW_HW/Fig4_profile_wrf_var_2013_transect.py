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

    ncfile  = Dataset(file_path, mode='r')

    Z       = getvar(ncfile, "z", timeidx=ALL_TIMES)
    Wa      = getvar(ncfile, "wa", units="m s-1", timeidx=ALL_TIMES)
    Ua      = getvar(ncfile, "ua", units="m s-1", timeidx=ALL_TIMES)
    T       = getvar(ncfile, 'th', units='degC', timeidx=ALL_TIMES)
    # T       = getvar(ncfile, 'temp', units='degC', timeidx=ALL_TIMES)
    S       = getvar(ncfile, 'QVAPOR', timeidx=ALL_TIMES)
    PBL     = getvar(ncfile, 'PBLH', timeidx=ALL_TIMES)

    print(PBL)

    return Z, Wa, Ua, T, S, PBL

def get_time_masked(Z,T,S,Ua,Wa,PBL,time_cood):

    z  = Z[time_cood,:,:,:]
    t  = T[time_cood,:,:,:]
    s  = S[time_cood,:,:,:]
    ua = Ua[time_cood,:,:,:]
    wa = Wa[time_cood,:,:,:]
    pbl= PBL[time_cood,:,:]

    return z, t, s, ua, wa, pbl

def get_vertcross(file_path, z, t, s, ua, wa, lat_slt, lon_min, lon_max, doy, seconds=None):

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.

    ntime       = np.shape(z)[0]

    ncfile      = Dataset(file_path)
    start_point = CoordPair(lat=lat_slt, lon=lon_min)
    end_point   = CoordPair(lat=lat_slt, lon=lon_max)

    # ****************** calc interpolation ******************
    t_out     = np.zeros((ntime, 36, 48))
    s_out     = np.zeros((ntime, 36, 48))
    ua_out    = np.zeros((ntime, 36, 48))
    wa_out    = np.zeros((ntime, 36, 48))

    for i in np.arange(ntime):

        print("i = ", i)

        t_crs  = vertcross(t[i,0:], z[i], wrfin=ncfile, start_point=start_point,
                            end_point=end_point, latlon=True, meta=True, autolevels=80)
        s_crs  = vertcross(s[i], z[i], wrfin=ncfile, start_point=start_point,
                            end_point=end_point, latlon=True, meta=True,autolevels=80)
        ua_crs = vertcross(ua[i], z[i], wrfin=ncfile, start_point=start_point,
                            end_point=end_point, latlon=True, meta=True,autolevels=80)
        wa_crs = vertcross(wa[i], z[i], wrfin=ncfile, start_point=start_point,
                            end_point=end_point, latlon=True, meta=True,autolevels=80)

        if i == 0:
            print("i == 0")
            loct = np.linspace(lon_min, lon_max, len(t_crs.coords['xy_loc']))
            print("loct", loct)
            vrt = np.arange(0, 3600., 100.)

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

def get_PBL(file_path, land_path, z, pbl, lat_slt, lon_min, lon_max):

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.

    ntime       = np.shape(z)[0]

    ncfile      = Dataset(file_path)
    start_point = CoordPair(lat=lat_slt, lon=lon_min)
    end_point   = CoordPair(lat=lat_slt, lon=lon_max)

    # ****************** calc interpolation ******************
    '''
    In vertcross, autolevels=100(default), then vertical profile is evenly spaced to 100 levels
    '''
    pbl_out     = np.zeros((ntime, 10, 48))
    pbl_4D      = np.expand_dims(pbl,axis=1).repeat(29,axis=1)

    landfile    = Dataset(land_path, mode='r')
    elev        = landfile.variables['Elevation_inst'][0]

    for i in np.arange(ntime):
        print(np.shape(z[i]))
        print(np.shape(pbl_4D[i]))

        pbl_crs = vertcross(pbl_4D[i]+elev, z[i], wrfin=ncfile, start_point=start_point,
                            end_point=end_point, latlon=True, meta=True, autolevels=10)
        print(np.shape(pbl_crs))
        print(pbl_crs)
        pbl_out[i] = pbl_crs
        pbl_crs    = None

    pbl_cross = np.nanmean(pbl_out[:,9,:], axis=0)
    print(pbl_cross)

    return pbl_cross

def get_WTD(file_path, land_path, time_s, time_e):

    ncfile      = Dataset(file_path)
    start_point = CoordPair(lat=lat_slt, lon=lon_min)
    end_point   = CoordPair(lat=lat_slt, lon=lon_max)
    z           = getvar(ncfile, "z")
    print(z)

    Time_s      = time_s - datetime(2000,1,1,0,0,0)
    Time_e      = time_e - datetime(2000,1,1,0,0,0)

    landfile    = Dataset(land_path, mode='r')
    Time        = nc.num2date(landfile.variables['time'][:],landfile.variables['time'].units,
                  only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time        = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)
    time_cood   = (time>=Time_s) & (time<Time_e)

    WTD         = landfile.variables['WaterTableD_tavg']
    wtd         = np.nanmean(WTD[time_cood,:,:],axis=0)
    wtd_3D      = np.expand_dims(wtd,axis=0).repeat(29,axis=0)
    print(np.shape(wtd_3D))
    print(np.shape(z))

    wtd_crs     = vertcross(wtd_3D, z, wrfin=ncfile, start_point=start_point,
                          end_point=end_point, latlon=True, meta=True, autolevels=10)
    print(wtd_crs)
    return wtd_crs[8:10,:]/1000.

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

def plot_profile_wrf_wind(file_paths, land_paths, time_s, time_e, message=None, lat_slt=36, lon_min=130, lon_max=160):

    # ****************** Get time coordiation ******************
    time_cood_all, time_cood_day, time_cood_night, doy_all, doy_day, doy_night = \
        get_time_cood(file_paths[0], time_s, time_e)
    print("time_cood_day",time_cood_day)
    print("time_cood_day",time_cood_night)

    # ****************** Get the WRF variables ******************
    Z1,Wa1,Ua1,T1,S1,PBL1 = read_wrf_var(file_paths[0])
    print("Z1", Z1)

    # ****************** Get time masked ******************
    z1_day, t1_day, s1_day, ua1_day, wa1_day, pbl1_day = \
                       get_time_masked(Z1, T1, S1, Ua1, Wa1, PBL1, time_cood_day)
    print("z1_day",z1_day)

    z1_night, t1_night, s1_night, ua1_night, wa1_night, pbl1_night = \
                       get_time_masked(Z1, T1, S1, Ua1, Wa1, PBL1, time_cood_night)

    # ****************** vertcross, interpolate and mean ******************
    seconds         = [6.*60.*60.,18.*60.*60.]
    t1_day_crs, s1_day_crs, ua1_day_crs, wa1_day_crs, xy_loc, vertical =\
        get_vertcross(file_paths[0], z1_day, t1_day, s1_day, ua1_day,
                      wa1_day, lat_slt, lon_min, lon_max, doy_day, seconds)

    pbl1_day_crs    = get_PBL(file_paths[0],  land_paths[0], z1_day, pbl1_day, lat_slt, lon_min, lon_max)

    seconds         = [18.*60.*60.,6.*60.*60.]
    t1_night_crs, s1_night_crs, ua1_night_crs, wa1_night_crs, xy_loc, vertical =\
        get_vertcross(file_paths[0], z1_night, t1_night, s1_night, ua1_night,
                      wa1_night, lat_slt, lon_min, lon_max, doy_night, seconds)
    pbl1_night_crs  = get_PBL(file_paths[0],  land_paths[0], z1_night, pbl1_night, lat_slt, lon_min, lon_max)

    if len(file_paths) > 1:
        # ****************** read second file ******************
        Z2,Wa2,Ua2,T2,S2,PBL2 = read_wrf_var(file_paths[1])

        z2_day, t2_day, s2_day, ua2_day, wa2_day, pbl2_day = \
                        get_time_masked(Z2,T2,S2,Ua2,Wa2,PBL2, time_cood_day)

        z2_night, t2_night, s2_night, ua2_night, wa2_night, pbl2_night = \
                        get_time_masked(Z2, T2, S2, Ua2, Wa2, PBL2, time_cood_night)

        seconds         = [6.*60.*60.,18.*60.*60.]
        t2_day_crs, s2_day_crs, ua2_day_crs, wa2_day_crs, xy_loc, vertical =\
            get_vertcross(file_paths[1], z2_day, t2_day, s2_day, ua2_day, wa2_day,
                          lat_slt, lon_min, lon_max, doy_day, seconds)

        pbl2_day_crs = get_PBL(file_paths[1], land_paths[1], z2_day, pbl2_day, lat_slt, lon_min, lon_max)

        seconds         = [18.*60.*60.,6.*60.*60.]
        t2_night_crs, s2_night_crs, ua2_night_crs, wa2_night_crs, xy_loc, vertical =\
            get_vertcross(file_paths[1], z2_night, t2_night, s2_night, ua2_night,
                          wa2_night, lat_slt, lon_min, lon_max, doy_night, seconds)

        pbl2_night_crs = get_PBL(file_paths[1], land_paths[1], z2_night, pbl2_night, lat_slt, lon_min, lon_max)

        # ****************** Get water table depth ******************
        wtd_crs   = get_WTD(file_paths[1], land_paths[1], time_s, time_e)


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

    # ===== set plot =====
    # color map
    color_map     = get_cmap("coolwarm")
    blue_map_neg  = truncate_colormap(color_map, minval=0., maxval=0.5)
    color_map     = get_cmap("coolwarm").reversed()
    blue_map_pos  = truncate_colormap(color_map, minval=0.5, maxval=1.)
    cmap          = color_map
    cmap1         = plt.cm.YlGnBu_r#get_cmap("Greens").reversed()

    # quiver scale
    if len(file_paths) > 1:
        scale = 1.
    else:
        scale = 20.

    # contour levels
    levels1   = [-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2]
    levels2   = [ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6]

    # Water table depth height
    wtd_hgt   = [0,300]

    # Day temperature
    contour   = ax[0,0].contourf(xy_loc, vertical, t_day_crs, levels=levels1, cmap=blue_map_neg, extend='both')
    cntr_wtd  = ax[0,0].contourf(xy_loc, wtd_hgt, wtd_crs, levels=np.arange(1,12,1), cmap=cmap1, extend='both')
    line1     = ax[0,0].plot(xy_loc,pbl1_day_crs,ls="-", color="black")
    line2     = ax[0,0].plot(xy_loc,pbl2_day_crs,ls="--", color="black")
    q         = ax[0,0].quiver(xy_loc[::3], vertical[::3], ua_day_crs[::3,::3],
                              wa_day_crs[::3,::3], angles='xy', scale_units='xy',
                              scale=scale, pivot='middle', color="white")
    ax[0,0].text(0.02, 0.95, "(a) Δθ$\mathregular{_{max}}$", transform=ax[0,0].transAxes, verticalalignment='top', bbox=props) # fontsize=14,
    ax[0,0].set_ylabel("Geopotential Height (m)")#, fontsize=12)

    # Night temperature
    contour   = ax[0,1].contourf(xy_loc, vertical, t_night_crs, levels=levels1, cmap=blue_map_neg, extend='both')
    cntr_wtd  = ax[0,1].contourf(xy_loc, wtd_hgt, wtd_crs, levels=np.arange(1,12,1), cmap=cmap1, extend='both')
    line1     = ax[0,1].plot(xy_loc,pbl1_night_crs,ls="-", color="black")
    line2     = ax[0,1].plot(xy_loc,pbl2_night_crs,ls="--", color="black")
    q         = ax[0,1].quiver(xy_loc[::3], vertical[::3], ua_night_crs[::3,::3],
                              wa_night_crs[::3,::3], angles='xy', scale_units='xy',
                              scale=scale, pivot='middle', color="white")
    ax[0,1].text(0.02, 0.95, "(b) Δθ$\mathregular{_{min}}$", transform=ax[0,1].transAxes, verticalalignment='top', bbox=props) # fontsize=14,
    cb_var    = fig.colorbar(contour, ax=ax[0], pad=0.01, orientation="vertical", aspect=20, shrink=0.88)
    cb_var.set_label('Δθ (${^o}$C)', loc='center') # rotation=270,


    # Day specific humidity
    color_map = get_cmap("coolwarm")
    cmap      = color_map.reversed()
    contour   = ax[1,0].contourf(xy_loc, vertical, s_day_crs*1000., levels=levels2, cmap=blue_map_pos, extend='both')
    cntr_wtd  = ax[1,0].contourf(xy_loc, wtd_hgt, wtd_crs, levels=np.arange(1,12,1), cmap=cmap1, extend='both')
    line1     = ax[1,0].plot(xy_loc,pbl1_day_crs,ls="-", color="black")
    line2     = ax[1,0].plot(xy_loc,pbl2_day_crs,ls="--", color="black")
    q         = ax[1,0].quiver(xy_loc[::3], vertical[::3], ua_day_crs[::3,::3],
                              wa_day_crs[::3,::3], angles='xy', scale_units='xy',
                              scale=scale, pivot='middle', color="white")
    ax[1,0].text(0.02, 0.95, "(c) Δq$\mathregular{_{day}}$", transform=ax[1,0].transAxes, verticalalignment='top', bbox=props) # fontsize=14,
    ax[1,0].set_xlabel("Longitude")#, fontsize=12)
    ax[1,0].set_ylabel("Geopotential Height (m)")#, fontsize=12)


    # Day specific humidity
    contour   = ax[1,1].contourf(xy_loc, vertical, s_night_crs*1000., levels=levels2, cmap=blue_map_pos, extend='both')
    cntr_wtd  = ax[1,1].contourf(xy_loc, wtd_hgt, wtd_crs, levels=np.arange(1,12,1), cmap=cmap1, extend='both')
    line1     = ax[1,1].plot(xy_loc,pbl1_night_crs,ls="-", color="black")
    line2     = ax[1,1].plot(xy_loc,pbl2_night_crs,ls="--", color="black")
    q         = ax[1,1].quiver(xy_loc[::3], vertical[::3], ua_night_crs[::3,::3],
                              wa_night_crs[::3,::3], angles='xy', scale_units='xy',
                              scale=scale, pivot='middle', color="white")
    
    ax[1,1].quiverkey(q,X=0.80, Y=2.1, U=scale, label=str(scale)+' m/s', labelpos='E', color="black")
    ax[1,1].text(0.02, 0.95, "(d) Δq$\mathregular{_{night}}$", transform=ax[1,1].transAxes, verticalalignment='top', bbox=props) # fontsize=14,
    ax[1,1].set_xlabel("Longitude")#, fontsize=12)

    cb_var    = fig.colorbar(contour, ax=ax[1], pad=0.01, orientation="vertical", aspect=20, shrink=0.88)
    cb_var.set_label('Δq (g kg$^{-1}$)', loc='center')

    # colorbar position
    position  = fig.add_axes([0.14, 0.04, 0.62, 0.02]) # [left, bottom, width, height]
    cb_wtd    = fig.colorbar(cntr_wtd, ax=ax, pad=0.07, cax=position, orientation="horizontal", aspect=40, shrink=0.8)
    
    cb_wtd.set_label('WTD (m)', loc='center',size=16)# rotation=270,
    cb_wtd.ax.tick_params(labelsize=12)

    fig.savefig("./plots/profile_wrf_Wind_"+message+"_transect", bbox_inches='tight', pad_inches=0.3)

if __name__ == "__main__":

    hw_name     = "hw2013_3Nov"
    lat_slt     = -26.
    lon_min     = 130.0
    lon_max     = 155.0

    # 2009
    if hw_name == "hw2009_3Nov":
        start_date= "20090122"
        end_date  = "20090213"
        # time_s = datetime(2009,1,28,5,0,0,0)
        # time_e = datetime(2009,1,28,6,59,0,0)
        time_s = datetime(2009,1,28,0,0,0,0)
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
    cpl_atmo_file     = '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/'+hw_name+'/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_'+start_date+'-'+end_date+'_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_'+start_date+'-'+end_date+'_fd'  # atmo output of wrf-cable run

    cpl_land_file     = '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/'+hw_name+'/ensemble_avg'
    cpl_land_file_gw  = cpl_land_file + '/LIS.CABLE.'+start_date+'-'+end_date+'_gw.nc' # land output of wrf-cable run
    cpl_land_file_fd  = cpl_land_file + '/LIS.CABLE.'+start_date+'-'+end_date+'_fd.nc' # land output of wrf-cable run

    file_paths        = [cpl_atmo_file_fd,cpl_atmo_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw
    land_paths        = [cpl_land_file_fd,cpl_land_file_gw]
    if len(file_paths) > 1:
        message = "GW-FD_"+str(time_s)+"-"+str(time_e)
    else:
        message = "GW_"+str(time_s)+"-"+str(time_e)

    plot_profile_wrf_wind(file_paths, land_paths, time_s, time_e, message=message, lat_slt=lat_slt, lon_min=lon_min, lon_max=lon_max)
