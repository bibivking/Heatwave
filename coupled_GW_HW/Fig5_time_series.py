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
from matplotlib.lines import Line2D
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
        file_tmp = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg/wrfout_20090122-20090213_gw"
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
        plt.savefig('./plots/check_pixels_rain_wtd_pft.png',dpi=300)

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

    '''
    output: time is AEST
    '''

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

def calc_diurnal_cycle(time,Var):

    seconds = [t.seconds for t in time]

    print("seconds",seconds)

    data = pd.DataFrame([t.seconds for t in time], columns=['seconds'])

    case_num = len(Var[0,:])

    for j in np.arange(case_num):
        data[str(j)] = Var[:,j]

    data_group = data.groupby(by=["seconds"]).mean()
    print(data_group)

    var = np.zeros([len(data_group),case_num])
    for j in np.arange(case_num):
        var[:,j] = data_group[str(j)].values
    print(var)

    return var

def calc_mean_min_max(var,ntime):

    var_mean      = np.zeros([2,ntime])
    var_min       = np.zeros([2,ntime])
    var_max       = np.zeros([2,ntime])

    var_mean[0,:] = np.nanmean(var[:,0:5],axis=1) # FD
    var_mean[1,:] = np.nanmean(var[:,5:11],axis=1) # GW
    print("np.shape(var[:,0:5])",np.shape(var[:,0:5]))

    var_min[0,:]  = np.nanmin(var[:,0:5],axis=1) # FD
    var_min[1,:]  = np.nanmin(var[:,5:11],axis=1) # GW

    var_max[0,:]  = np.nanmax(var[:,0:5],axis=1) # FD
    var_max[1,:]  = np.nanmax(var[:,5:11],axis=1) # GW

    return var_mean,var_min,var_max

def plot_time_series_3hws( path, case_names, periods, time_ss, time_es, seconds=None,loc_lats=None, loc_lons=None,
                      message=None, rain_val=None, wtd_val=None, pft_val=None, is_diurnal=True, is_envelop=True):

    print("======== In plot_time_series =========")

    # ============== Set the plot ==============
    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=[12,12], sharex=True, squeeze=True) # sharex=True,sharey='row',
    plt.subplots_adjust(left=0.06,top=0.98,right=0.9, wspace=0.155, hspace=0)

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

    # set the box type of sequence number
    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')

    # ======================= Set colormap =======================
    cmap      = plt.cm.seismic

    colors    = [ ["orange","green"],  # pre-heatwave
                    ["brown","blue"]     ] # heatwave
    labels    = [ ["pre-hw FD","pre-hw GW"],  # pre-heatwave
                    ["hw FD","hw GW"]         ] # heatwave
    orders    = ['(a)','(b)','(c)','(d)','(e)','(f)']

    label_x   = ["T$\mathregular{_{2m}}$ ($^{o}$C)","ABL (m)"]
    label_y   = ["2009","2013","2019"]
    loc_y     = [0.55,0.5,0.45]

    # ============== read heatwave information ==============
    for hw_num in np.arange(3):
        case_name = case_names[hw_num]
        period    = periods[hw_num]

        # ============== set up file names ==============
        if case_name == "hw2009_3Nov":
            rst_dates = ["20090117","20090118","20090119","20090120","20090121" ]
        elif case_name == "hw2013_3Nov":
            rst_dates = ["20121224","20121225","20121226","20121227","20121228" ]
        elif case_name == "hw2019_3Nov":
            rst_dates = ["20190103","20190104","20190105","20190106","20190107" ]

        file_path  = path + case_name
        file_paths = None
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
        case_sum  = len(file_paths)

        # ============== loop different location shallow WTD vs deep WTD ==============
        loc_num      = 0

        loc_lat      = loc_lats[loc_num]
        loc_lon      = loc_lons[loc_num]
        print("lat=",loc_lat,"lon=",loc_lon)

        x_ticks      = []
        x_ticklabels = []
        nperiods     = len(time_ss[hw_num][:])

        # ============== loop different periods (pre/during heatwaves) ==============
        for cyc_num in np.arange(nperiods):
            time_s    = time_ss[hw_num][cyc_num]
            time_e    = time_es[hw_num][cyc_num]
            print("period=",time_s,"-",time_e)

            # ============== read variables ==============
            time, tmp = read_ensembles(file_paths[0], "PBL", time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon)
            tmp       = None

            ntime     = len(time) # number of total selected timestep
            tair      = np.zeros([ntime,case_sum])
            pbl       = np.zeros([ntime,case_sum])
            # hc        = np.zeros([ntime,case_sum])

            # ============== get mask =============
            is_mask     = (rain_val is not None) | (wtd_val is not None) | (pft_val is not None)
            check_pixel = False

            if is_mask:
                land_file  = path + case_name + "/ensemble_avg/LIS.CABLE."+period+"_gw.nc"
                print("land_file=",land_file)
                mask_map = get_mask(land_file, time_s, time_e, rain_val, wtd_val, pft_val, check_pixel=check_pixel)

            # ============== read variables =============
            for case_num in np.arange(case_sum):
                file_path = file_paths[case_num]
                time, tair[:,case_num]= read_ensembles(file_path, "Tair", time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon)
                time, pbl[:,case_num] = read_ensembles(file_path, "PBL", time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon)
                # hc_tmp[:,case_num]    = read_heat_content(file_path, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon)

            # ============== calc diurnal ==============
            if is_diurnal:
                tair = calc_diurnal_cycle(time,tair)
                pbl  = calc_diurnal_cycle(time,pbl)
                # hc   = calc_diurnal_cycle(time,hc)

            # ============== plotting ==============
            if cyc_num%2 ==0:
                # per heatwave
                nt  = len(tair[:,0])
                x   = np.arange(nt) + (nt+1)*int(cyc_num/2)

            if is_envelop:

                tair_mean,tair_min,tair_max = calc_mean_min_max(tair,nt)
                pbl_mean,pbl_min,pbl_max    = calc_mean_min_max(pbl,nt)
                # hc_mean,hc_min,hc_max       = calc_mean_min_max(hc,nt)

                for fd_gw in np.arange(2):

                    color = colors[cyc_num%2][fd_gw]
                    label = labels[cyc_num%2][fd_gw]

                    # Tair
                    var   = tair_mean[fd_gw,:]
                    x_max = np.argmax(var)
                    ax[hw_num,0].plot(x, var, color=color, lw=1.0, label=label)
                    ax[hw_num,0].fill_between(x, tair_min[fd_gw,:], tair_max[fd_gw,:], alpha=0.5, facecolor=color)
                    ax[hw_num,0].axvline(x=x[x_max], color=color,lw=1.0, alpha=0.4, linestyle='--')

                    # PBL
                    var   = pbl_mean[fd_gw,:]
                    x_max = np.argmax(var)
                    ax[hw_num,1].plot(x, var, color=color, lw=1.0, label=label)
                    ax[hw_num,1].fill_between(x, pbl_min[fd_gw,:], pbl_max[fd_gw,:], alpha=0.5, facecolor=color)
                    # ax[1].axvline(x=x[x_max], color=color,lw=1.0, alpha=0.4, linestyle='--')

            if cyc_num%2 == 0:
                if is_diurnal:
                    x_ticks      = np.concatenate((x_ticks, x[::6]), axis=0)
                    x_ticklabels = np.concatenate((x_ticklabels, ['0:00','6:00','12:00','18:00']), axis=0)
                else:
                    # ax.set_xlim([0,nt])
                    x_ticks      = np.arange(0,max(x),12)
                    x_ticklabels = ['-1d 0:00','-1d 12:00','1d 0:00','1d 12:00','2d 0:00','2d 12:00',
                                    '3d 0:00','3d 12:00','4d 0:00','4d 12:00','+1d 0:00','+1d 12:00']

        print(x_ticks)
        print(x_ticklabels)

        # =============== add label & order ===============

        # set x label/title
        if hw_num == 0:
            ax[0,0].set_title(label_x[0],fontsize=16)
            ax[0,1].set_title(label_x[1],fontsize=16)

        # set y label
        ax[hw_num,0].text(-0.1, loc_y[hw_num], label_y[hw_num], va='bottom', ha='center',
                            rotation='vertical', rotation_mode='anchor',
                            transform=ax[hw_num,0].transAxes)

        ax[hw_num,0].set_xticks(x_ticks)
        ax[hw_num,0].set_xticklabels(x_ticklabels,fontdict={'fontsize':14})
        ax[hw_num,0].tick_params(axis='y', labelsize=14)
        ax[hw_num,0].text(0.02, 0.95, orders[hw_num*2], transform=ax[hw_num,0].transAxes, verticalalignment='top', bbox=props,fontdict={'fontsize':18})

        ax[hw_num,1].set_xticks(x_ticks)
        ax[hw_num,1].set_xticklabels(x_ticklabels,fontdict={'fontsize':14})
        ax[hw_num,1].tick_params(axis='y', labelsize=14)
        ax[hw_num,1].text(0.02, 0.95, orders[hw_num*2+1], transform=ax[hw_num,1].transAxes, verticalalignment='top', bbox=props,fontdict={'fontsize':18})

        # =============== add legend ===============
        custom_lines = [Line2D([0], [0], color="orange", lw=1),
                        Line2D([0], [0], color="green", lw=1),
                        Line2D([0], [0], color="brown", lw=1),
                        Line2D([0], [0], color="blue", lw=1)]

        fig.legend(custom_lines, ['pre-hw FD', 'pre-hw GW', 'hw FD', 'hw GW'],
                loc='upper right', bbox_to_anchor=(0.84, 0.95), frameon=False, fontsize=14)
        # fig.tight_layout()

        # ============ savefig ============
        if is_diurnal:
            message = 'durinal_cycle_Tmax_PBL_' +case_name
        else:
            message = 'time_series_Tmax_PBL_' +case_name

        # if loc_lat != None:
        #     message = message + "_lat="+str(loc_lat) + "_lon="+str(loc_lon)

        # fig.tight_layout()
    plt.savefig('./plots/'+message+'.png',dpi=300,bbox_inches='tight', pad_inches=0.3)


if __name__ == "__main__":

    # #######################
    #        Set decks      #
    # #######################

    path       = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"

    case_names = [
                  "hw2009_3Nov",
                  "hw2013_3Nov",
                  "hw2019_3Nov"]

    periods    = ["20090122-20090213",
                  "20121229-20130122",
                  "20190108-20190130"]

    #              pre-heatwave                    heatwave
    time_ss      = [ [datetime(2009,1,24,0,0,0),   datetime(2009,1,28,0,0,0),
                      datetime(2009,2,1,0,0,0),    datetime(2009,2,5,0,0,0)   ],
                     [datetime(2013,1,1,0,0,0),    datetime(2013,1,5,0,0,0),
                      datetime(2013,1,15,0,0,0),   datetime(2013,1,17,0,0,0)    ],
                     [datetime(2019,1,9,0,0,0),    datetime(2019,1,14,0,0,0),
                      datetime(2019,1,20,0,0,0),   datetime(2019,1,22,0,0,0)   ] ]

    time_es      = [ [datetime(2009,1,27,23,59,0), datetime(2009,1,31,23,59,0),
                      datetime(2009,2,4,23,59,0),  datetime(2009,2,8,23,59,0)  ],
                     [datetime(2013,1,4,23,59,0),  datetime(2013,1,8,23,59,0),
                      datetime(2013,1,16,23,59,0), datetime(2013,1,18,23,59,0)  ],
                     [datetime(2019,1,13,23,59,0), datetime(2019,1,18,23,59,0),
                      datetime(2019,1,21,23,59,0), datetime(2019,1,26,23,59,0) ] ]


    loc_lats   = [[-36.5,-35.5],
                  [-36.5,-35.5]]


    loc_lons   = [[148.,149.],
                  [141.5,142.5]]


    # loc_lons   = [[148.5,149.5],
    #               [141.5,142.5]]

    plot_time_series_3hws(path, case_names, periods, time_ss, time_es,
                loc_lats=loc_lats, loc_lons=loc_lons, is_diurnal=True, is_envelop=True)
