#!/usr/bin/python

'''
Functions:
1. plot multi lines in one figure
'''

from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from matplotlib import cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.cm import get_cmap
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_scale_offline
from common_utils import *

def get_mask(land_file, time_s, time_e, rain_val=None, wtd_val=None, pft_val=None, check_pixel=False):

    # get pixels wants to keep in the whole domain --- Note that this function weren't tested

    if rain_val is not None:
        Time, rain= read_var(land_file, "Rainf_tavg", lat_name="lat", lon_name="lon")
        rain_mean = spital_var(UTC_to_AEST(Time), rain, time_s, time_e)
        rain_mask = np.ones(np.shape(rain_mean), dtype=bool)
        rain_mask = np.where( rain_mean*24.*3600. < rain_val, True, False)
        print("(rain_mask == True).sum()")
        print((rain_mask == True).sum())
        mask_map = rain_mask

    if wtd_val is not None:
        Time, wtd = read_var(land_file, "WaterTableD_tavg", lat_name="lat", lon_name="lon")
        wtd_mean = spital_var(UTC_to_AEST(Time), wtd, time_s, time_e)
        wtd_mask = np.ones(np.shape(wtd_mean), dtype=bool)
        wtd_mask = np.where(np.all([(wtd_mean/1000.) >=wtd_val[0], (wtd_mean/1000.) <wtd_val[1]],axis=0), True, False)
        print("(wtd_mask == True).sum()")
        print((wtd_mask == True).sum())
        if mask_map is not None:
            mask_map = np.all([mask_map,wtd_mask],axis=0)
        else:
            mask_map = wtd_mask

    if pft_val is not None:
        pft      = tree_mask(land_file,"Landcover_inst")
        pft_mask = np.ones(np.shape(pft), dtype=bool)
        if len(pft_val) == 1:
            pft_mask = np.where( pft == pft_val[0], True, False)
        elif (len(pft_val) == 2) :
            pft_mask = np.where( np.all([pft >=pft_val[0], pft<=pft_val[1]],axis=0), True, False)
        if mask_map is not None:
            mask_map = np.all([mask_map,pft_mask],axis=0)
        else:
            mask_map = pft_mask

    # check pixel cover
    if check_pixel:

        # check mask
        file_tmp = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg/wrfout_20090122-20090213_gw"
        file = nc.Dataset(file_tmp, mode='r')
        p1   = getvar(file, "pressure")
        lon  = file.variables['XLONG'][0,:,:]
        lat  = file.variables['XLAT'][0,:,:]

        fig1 = plt.figure(figsize=(12,9))
        # Get the lat/lon coordinates

        ax1 = plt.axes(projection=ccrs.PlateCarree())
        ax1.set_extent([110,155,-45,-10])
        ax1.coastlines(resolution="50m",linewidth=1)

        # Add gridlines
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top  = False
        gl.ylabels_right= False
        gl.xlines       = True
        gl.xlocator     = mticker.FixedLocator(np.arange(110,155,1))
        gl.ylocator     = mticker.FixedLocator(np.arange(-45,-10,1))
        gl.xformatter   = LONGITUDE_FORMATTER
        gl.yformatter   = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black','rotation': 90}
        gl.ylabel_style = {'size':10, 'color':'black'}

        mask = np.where(np.isnan(mask_map), 0, 1)
        plt.contourf(lon, lat, mask, levels=[-0.1,0.1,1.1], transform=ccrs.PlateCarree(), cmap=plt.cm.jet)
        plt.savefig('./plots/figures/check_pixels_rain_wtd_pft.png',dpi=300)

    return mask_map

def air_density(T,QV,P):

    # source: https://mailman.ucar.edu/pipermail/ncl-talk/2019-August/015431.html

    '''
    Ideal gas law:  p = (density)*RD*(temperature)
                        Hence:
                        density = p/(RD*temperature)
    p   => Pa
    qv  => kg/kg  :QVAPOR
    qc  => kg/kg  :QCLOUD
    t   => K
    '''

    RD      = 287.04         # [J/(K-kg)] Gas Constant dry air
    TV      = T*(1+0.609*QV)
    density = P/(RD*TV)

    return density

def calc_atmospheric_heat_content(t,density):

    """
    Quantify ocean heat content upto certain depths.
    param T: interpolation of potential temperature on the vector grid (v grid)    [z, jj, ji]
    param mask: land-sea mask on vector grid               [z, jj, ji]
    param e1t: width of the grid box in zonal direction on T grid    [jj, ji]
    param e2t: height of the grid box in zonal direction on T grid   [jj, ji]
    param e3t_0: depth of each grid box                    [z]
    param e3t_adjust: adjustment of grid box depth based on the definition of partial grid  [z, jj, ji]

    return: arrays of OHC upto certain depth
    rtype: numpy arrays
    """

    # define the constant:
    constant = {
                'g' : 9.80616,      # gravititional acceleration [m / s2]
                'cp': 1006,         # heat capacity of air [J/(Kg*C)]
                'rho': 1027,        # sea water density [Kg/m3]
                }

    ntime = len(t[:,0])

    # calculate heat flux at each grid point
    AHC = np.zeros(np.shape(t))

    for i in np.arange(len(t[0,:])):
        AHC[:,i] = density[:,i]*constant['cp']*t[:,i]*20.
        AHC_sum  = np.sum(AHC,axis=1)/1e+12
    print(AHC_sum)
    
    return AHC_sum

def read_heat_content(file_path, time_s, time_e, loc_lat=None, loc_lon=None, mask_map=None):

    """
    Quantify atmospheric heat content upto certain depths.
    param T: interpolation of potential temperature on the vector grid (v grid)    [z, jj, ji]
    param e1t: width of the grid box in zonal direction on T grid    [jj, ji]
    param e2t: height of the grid box in zonal direction on T grid   [jj, ji]
    param e3t_0: depth of each grid box                    [z]
    param e3t_adjust: adjustment of grid box depth based on the definition of partial grid  [z, jj, ji]

    return: arrays of OHC upto certain depth
    rtype: numpy arrays
    """

    time   = read_wrf_time(file_path)
    ntime  = len(time)
    print(time)
    
    var    = read_wrf_hgt_var(file_path, "temp", "K", height=10, loc_lat=loc_lat, loc_lon=loc_lon, p_hgt="hgt")
    nlat   = len(var[0,:,0])
    nlon   = len(var[0,0,:])
    var    = None
    
    # PBLH   = []
    T      = np.zeros([ntime,5,nlat,nlon])
    QV     = np.zeros([ntime,5,nlat,nlon])
    P      = np.zeros([ntime,5,nlat,nlon])
    cnt    = 0
    
    for height in np.arange(10, 110,20):
        # PBLH.append(read_wrf_hgt_var(file_path, "PBLH", "m", height=height, loc_lat=loc_lat, loc_lon=loc_lon, p_hgt="hgt")) 
        T[:,cnt,:,:]  = read_wrf_hgt_var(file_path, "temp", "K", height=height, loc_lat=loc_lat, loc_lon=loc_lon, p_hgt="hgt")
        QV[:,cnt,:,:] = read_wrf_hgt_var(file_path, "QVAPOR", height=height, loc_lat=loc_lat, loc_lon=loc_lon, p_hgt="hgt") # kg kg-1
        P[:,cnt,:,:]  = read_wrf_hgt_var(file_path, "p", "Pa",height=height, loc_lat=loc_lat, loc_lon=loc_lon, p_hgt="hgt")
        cnt = cnt + 1 
    print(np.shape(T))
    
    # pblh  =  time_series_var(time,PBLH,time_s,time_e)
    t     =  time_series_var(time,T,time_s,time_e)
    qv    =  time_series_var(time,QV,time_s,time_e)
    p     =  time_series_var(time,P,time_s,time_e)

    density = air_density(t,qv,p)
    hc      = calc_atmospheric_heat_content(t,density)
    print(hc)
    
    return hc

def read_ensembles(file_path, var_name, time_s, time_e, loc_lat=None, loc_lon=None, mask_map=None):

    print("In read_ensembles")

    Time    = read_wrf_time(file_path)
    print(var_name)
    if var_name == "Tair":
        Var      = read_wrf_surf_var(file_path, "T2", loc_lat=loc_lat, loc_lon=loc_lon, mask_map=mask_map) - 273.15
        time,var = time_series_var(Time,Var,time_s,time_e)
    elif var_name == "PBL":
        Var       = read_wrf_surf_var(file_path, "PBLH", loc_lat=loc_lat, loc_lon=loc_lon, mask_map=mask_map)
        time,var  = time_series_var(Time,Var,time_s,time_e)

    print(var)
    return time, var

def read_heat_thrhld( time_s,time_e, loc_lat=None, loc_lon=None):
    filename = "./nc_file/AWAP_tmean_1970_2019_90-th.nc"

    # file = nc.Dataset(filename, mode='r')
    time, tmean_90 = read_var(filename, "tmean", loc_lat, loc_lon, "latitude", "longitude")
    time_cood      = time_mask(time, time_s, time_e)
    tmean          = np.nanmean(tmean_90[time_cood])

    return tmean

def plot_time_series( path, case_names, periods, time_ss, time_es, seconds=None,loc_lats=None, loc_lons=None,
                      message=None, rain_val=None, wtd_val=None, pft_val=None):

    print("======== In plot_time_series =========")

    # ============== Set the plot ==============
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=[8,4],sharey=True, squeeze=True) # sharex=True, 
    plt.subplots_adjust(wspace=0.0, hspace=0.0)
    
    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
    print(ax)

    # ============== read heatwave information ==============
    hw_num    = 0
    case_name = case_names[hw_num]
    period    = periods[hw_num]
    time_s    = time_ss[hw_num]
    time_e    = time_es[hw_num]
    loc_lat   = loc_lats[hw_num]
    loc_lon   = loc_lons[hw_num]

    if case_name == "hw2009_3Nov":
        rst_dates = ["20090117","20090118","20090119","20090120","20090121" ]
    elif case_name == "hw2013_3Nov":
        rst_dates = ["20121224","20121225","20121226","20121227","20121228" ]
    elif case_name == "hw2019_3Nov":
        rst_dates = ["20190103","20190104","20190105","20190106","20190107" ]

    # ============== set up file names ==============
    file_path  = path + case_name
    file_paths = [ file_path + "/fd_rst_" + rst_dates[0] + "/WRF_output/wrfout_" + period,
                   file_path + "/fd_rst_" + rst_dates[1] + "/WRF_output/wrfout_" + period,
                   file_path + "/fd_rst_" + rst_dates[2] + "/WRF_output/wrfout_" + period,
                   file_path + "/fd_rst_" + rst_dates[3] + "/WRF_output/wrfout_" + period,
                   file_path + "/fd_rst_" + rst_dates[4] + "/WRF_output/wrfout_" + period,
                   file_path + "/gw_rst_" + rst_dates[0] + "/WRF_output/wrfout_" + period,
                   file_path + "/gw_rst_" + rst_dates[1] + "/WRF_output/wrfout_" + period,
                   file_path + "/gw_rst_" + rst_dates[2] + "/WRF_output/wrfout_" + period,
                   file_path + "/gw_rst_" + rst_dates[3] + "/WRF_output/wrfout_" + period,
                   file_path + "/gw_rst_" + rst_dates[4] + "/WRF_output/wrfout_" + period ]

    # ============== read variables ==============
    case_sum = len(file_paths)
    
    time, tmp = read_ensembles(file_paths[0], "PBL", time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon)
    tmp      = None
    ntime    = len(time)
    pbl      = np.zeros([ntime,case_sum])
    hc       = np.zeros([ntime,case_sum])
    tair     = np.zeros([ntime,case_sum])


    # ============== get mask =============
    is_mask     = (rain_val is not None) | (wtd_val is not None) | (pft_val is not None)
    check_pixel = False

    if is_mask:
        land_file  = path + case_name + "/ensemble_avg/LIS.CABLE."+period+"_gw.nc"
        print("land_file=",land_file)
        mask_map = get_mask(land_file, time_s, time_e, rain_val, wtd_val, pft_val, check_pixel=check_pixel)

    for case_num in np.arange(case_sum):
        file_path = file_paths[case_num]
        time, tair[:,case_num]= read_ensembles(file_path, "Tair", time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon)
        time, pbl[:,case_num] = read_ensembles(file_path, "PBL", time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon)
        hc_tmp[:,case_num] = read_heat_content(file_path, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon)
        
    print("time = ", time )
    print("np.shape(tair) = ", np.shape(tair) )
    print("np.shape(pbl) = ", np.shape(pbl) )
    print("np.shape(hc) = ", np.shape(hc) )


    # ============== plotting ==============
    labels = ["fd-10d","fd-9d","fd-8d","fd-7d","fd-6d","gw-10d","gw-9d","gw-8d","gw-7d","gw-6d"]
    alpha  = [1.0,0.9,0.8,0.7,0.6,1.0,0.9,0.8,0.7,0.6]
    # colors = []
    for case_num in np.arange(case_sum):
        if case_num < 5:
            plot1 = ax[0].plot(tair[:,case_num],  color="darkred", alpha = alpha[case_num], label=labels[case_num]) # colors="Reds",
            plot2 = ax[1].plot(pbl[:,case_num], color="darkred", alpha = alpha[case_num], label=labels[case_num]) # colors="Reds",
            plot3 = ax[2].plot(hc[:,case_num], color="darkred", alpha = alpha[case_num], label=labels[case_num])
        else:
            plot1 = ax[0].plot(tair[:,case_num],  color="darkblue", alpha = alpha[case_num], label=labels[case_num]) 
            plot2 = ax[1].plot(pbl[:,case_num], color="darkblue", alpha = alpha[case_num], label=labels[case_num]) # colors="Reds",
            plot3 = ax[2].plot(hc[:,case_num], color="darkred", alpha = alpha[case_num], label=labels[case_num])
    
    # heat_thrhld = read_heat_thrhld( time_s,time_e, loc_lat=loc_lat, loc_lon=loc_lon)
    # print(heat_thrhld)
    # ax[0].axhline(y=heat_thrhld, color="gray", linestyle='--')
    ax[0].text(0.02, 0.95, "(a) T", transform=ax[0].transAxes, verticalalignment='top', bbox=props)
    ax[1].text(0.02, 0.95, "(b) PBL", transform=ax[1].transAxes, verticalalignment='top', bbox=props)
    ax[2].text(0.02, 0.95, "(c) HC", transform=ax[1].transAxes, verticalalignment='top', bbox=props)
    
    x_ticks      = np.arange(ntime,12)
    print(ntime)
    x_ticklabels = ['-1d 0:00','-1d 12:00','1d 0:00','1d 12:00','2d 0:00','2d 12:00',
                    '3d 0:00','3d 12:00','4d 0:00','4d 12:00','+1d 0:00','+1d 12:00']
    ax[2].set_xticks(x_ticks)
    ax[2].set_xticklabels(x_ticklabels)
    # plt.legend(plot1)
    # ax.legend()

    
    # ax.set_xlim([np.min(var1*scale,var2*scale), np.max(var1*scale,var2*scale)])
    # ax.plot(t1, var*scale, alpha=0.5)
    # ax.set_ylabel('mm')
    # ax.set_title(var_name)

    # fig.tight_layout()
    # if message == None:
    #     message = var_name
    # else:
    #     message = message + "_" + var_name
    # if loc_lat != None:
    message = message + "_lat="+str(loc_lat) + "_lon="+str(loc_lon)

    plt.savefig('./plots/figures/test_time_series_Tmax_PBL_'+message+'.png',dpi=300)



if __name__ == "__main__":

    # #######################
    #        Set decks      #
    # #######################

    path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"

    case_names = ["hw2009_3Nov",
                  "hw2013_3Nov",
                  "hw2019_3Nov"]

    periods    = ["20090122-20090213",
                  "20121229-20130122",
                  "20190108-20190130"]

    time_ss    = [datetime(2009,1,28,0,0,0),
                  datetime(2013,1,4,0,0,0),
                  datetime(2019,1,14,0,0,0)]

    time_es    = [datetime(2009,1,28,23,59,0),
                  datetime(2013,1,22,23,59,0),
                  datetime(2019,1,30,23,59,0)]

    seconds    = [6.*60.*60.,18.*60.*60.]

    loc_lats   = [[-37,-35],
                  [-29,-28]]

    loc_lons   = [[148.,150.],
                  [138.,139.]]

    if seconds == None:
        time_spell = "all"
    elif seconds[0] < seconds[1]:
        time_spell = "day"
    elif seconds[0] > seconds[1]:
        time_spell = "night"

    message   = time_spell

    plot_time_series(path, case_names, periods, time_ss, time_es, seconds,
                loc_lats=loc_lats, loc_lons=loc_lons, message=message )
