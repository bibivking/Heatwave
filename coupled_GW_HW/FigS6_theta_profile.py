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
from common_utils import *

def get_day_night_cood(file_path, time_s ,time_e):

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

    # for daytime
    seconds       = [12.*60.*60.,
                     16.*60.*60.]

    time_cood_day = []
    for j in np.arange(len(time)):
        if_day    = (time[j].seconds >= seconds[0]) & (time[j].seconds < seconds[1])
        time_cood_day.append((time[j]>=Time_s) & (time[j]<Time_e) & if_day)

    # for nighttime
    seconds       = [0.*60.*60.,
                     4.*60.*60.]

    time_cood_night = []
    for j in np.arange(len(time)):
        if_night  = (time[j].seconds >= seconds[0]) & (time[j].seconds < seconds[1])
        time_cood_night.append((time[j]>=Time_s) & (time[j]<Time_e) & if_night)

    return time_cood_day, time_cood_night

def read_hgt_theta(file_path, time_s, time_e, loc_lat=None, loc_lon=None):

    # read var
    ncfile  = Dataset(file_path, mode='r')
    Z       = getvar(ncfile, "height_agl", timeidx=ALL_TIMES)
    T       = getvar(ncfile, 'th', units='degC', timeidx=ALL_TIMES)

    # get time time
    time_cood_day, time_cood_night = get_day_night_cood(file_path, time_s, time_e)
    # mask by time
    z_day    = Z[time_cood_day,:,:,:]
    z_night  = Z[time_cood_night,:,:,:]
    t_day    = T[time_cood_day,:,:,:]
    t_night  = T[time_cood_night,:,:,:]

    print(np.shape(z_day))
    print(np.shape(z_night))
    print(np.shape(t_day))
    print(np.shape(t_night))

    # get spitial mask
    mask   = mask_by_lat_lon(file_path, loc_lat, loc_lon, 'XLAT', 'XLONG')
    # mask by location
    for time in np.arange(len(z_day[:,0,0,0])):
        for lvl in np.arange(len(z_day[0,:,0,0])):
            z_day[time,lvl,:,:]    = np.where(mask,z_day[time,lvl,:,:],np.nan)
            z_night[time,lvl,:,:]  = np.where(mask,z_night[time,lvl,:,:],np.nan)
            t_day[time,lvl,:,:]    = np.where(mask,t_day[time,lvl,:,:],np.nan)
            t_night[time,lvl,:,:]  = np.where(mask,t_night[time,lvl,:,:],np.nan)

    print(t_day[0,:,:,:])
    hgt_day    = np.nanmean(z_day,axis=(0,2,3))
    hgt_night  = np.nanmean(z_night,axis=(0,2,3))
    theta_day  = np.nanmean(t_day,axis=(0,2,3))
    theta_night= np.nanmean(t_night,axis=(0,2,3))
    print("hgt_day, hgt_night, theta_day, theta_night",hgt_day, hgt_night, theta_day, theta_night)
    return hgt_day, hgt_night, theta_day, theta_night

def plot_theta_profile( path, case_names, periods, time_ss, time_es, seconds=None, loc_lats=None, loc_lons=None, message=None):

    print("======== In plot_theta_profile =========")

    # ============== read heatwave information ==============
    # hw_num    = 1 # 2009

    for hw_num in np.arange(3):
        case_name = case_names[hw_num]
        period    = periods[hw_num]

        # ============== set up file names ==============
        cpl_atmo_file     = path+case_name+'/ensemble_avg'
        cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_'+periods[hw_num]+'_gw'  # atmo output of wrf-cable run
        cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_'+periods[hw_num]+'_fd'  # atmo output of wrf-cable run
        file_paths        = [cpl_atmo_file_fd,cpl_atmo_file_gw]
        case_sum          = len(file_paths)

        # ================ Start to process data ====================
        nperiods          = len(time_ss[hw_num][:])
        nlvl              = 29

        # ============== Set the plot ==============
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[15,4.5], sharex=True, squeeze=True) # sharex=True,sharey='row',
        plt.subplots_adjust(left=0.06,top=0.98,right=0.9, wspace=0.14, hspace=0)

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

        for cyc_num in np.arange(nperiods):

            # ============== loop different periods (pre/during heatwaves) ==============
            time_s    = time_ss[hw_num][cyc_num]
            time_e    = time_es[hw_num][cyc_num]
            print("period=",time_s,"-",time_e)

            pre_HW    = cyc_num%2      # 0 is pre-hw, 1 is hw
            HW_evt    = int(cyc_num/2) # 0 is hw1, 1 is hw2

            for loc_num in np.arange(1):

                hgt   = np.zeros([2,case_sum,nlvl]) # day/night x FW/GW x atmospheric layers
                theta = np.zeros([2,case_sum,nlvl]) # day/night x FW/GW x atmospheric layers

                # ============== loop different location shallow WTD vs deep WTD ==============
                loc_lat      = loc_lats[loc_num]
                loc_lon      = loc_lons[loc_num]
                print("lat=",loc_lat,"lon=",loc_lon)

                # ============== read variables =============
                for case_num in np.arange(case_sum):
                    file_path = file_paths[case_num]

                    hgt[0,case_num,:], hgt[1,case_num,:], theta[0,case_num,:], theta[1,case_num,:] = \
                        read_hgt_theta(file_path, time_s, time_e, loc_lat=loc_lats[loc_num], loc_lon=loc_lons[loc_num])

                # if loc_num == 1:
                #     theta[:,:,:] = theta[:,:,:] + 60

                if pre_HW == 0:
                    # pre_HW day
                    label     = "pre-hw"
                    # day
                    plot  = ax[HW_evt].plot(theta[0,0,:], hgt[0,0,:], color="orange", linestyle='-', lw=1.0, label=label+" FD")
                    plot  = ax[HW_evt].plot(theta[0,1,:], hgt[0,1,:], color="green", linestyle='-', lw=1.0, label=label+" GW")
                    # night
                    plot  = ax[HW_evt].plot(theta[1,0,:]+ 25., hgt[1,0,:], color="orange", linestyle='-', lw=1.0, label=label+" FD")
                    plot  = ax[HW_evt].plot(theta[1,1,:]+ 25., hgt[1,1,:], color="green", linestyle='-', lw=1.0, label=label+" GW")
                else:
                    # HW day
                    label     = "hw"
                    # day
                    plot  = ax[HW_evt].plot(theta[0,0,:], hgt[0,0,:], color="red", linestyle='-', lw=1.0, label=label+" FD")
                    plot  = ax[HW_evt].plot(theta[0,1,:], hgt[0,1,:], color="blue", linestyle='-', lw=1.0, label=label+" GW")
                    # night
                    plot  = ax[HW_evt].plot(theta[1,0,:]+ 25., hgt[1,0,:], color="red", linestyle='-', lw=1.0, label=label+" FD")
                    plot  = ax[HW_evt].plot(theta[1,1,:]+ 25., hgt[1,1,:], color="blue", linestyle='-', lw=1.0, label=label+" GW")

        ax[0].axvline(x=45, color=almost_black,lw=1.0, linestyle='-')
        # ax[0].axvline(x=65, color="gray", linestyle='-')
        # ax[0].axvline(x=85, color="gray", linestyle='--')

        ax[1].axvline(x=45, color=almost_black,lw=1.0, linestyle='-')
        # ax[1].axvline(x=65, color="gray", linestyle='-')
        # ax[1].axvline(x=85, color="gray", linestyle='--')

        # =============== add coordinates and order ===============
        x_ticks      = [25,30,35,40, # loc 1 day
                        50,55,60,65 # loc 1 night
                        # 70,75,80, # loc 2 day
                        # 90,95,100  # loc 2 day
                        ]
        x_ticklabels = ['25','30','35','40',  # loc 1 day
                        '25','30','35','40',  # loc 1 night
                        # '30','35','40',  # loc 2 day
                        # '30','35','40',  # loc 2 day
                        ]
        print(x_ticks)
        print(x_ticklabels)
        ax[0].set_xticks(x_ticks)
        ax[0].set_xticklabels(x_ticklabels)
        ax[0].set_xlim(20,70)
        ax[0].set_ylim(0,2500)
        ax[0].text(0.02, 0.95, "(a)", transform=ax[0].transAxes, verticalalignment='top', bbox=props)
        ax[0].text(0.02, 0.08, "afternoon", transform=ax[0].transAxes, verticalalignment='top', bbox=props)
        ax[0].text(0.9 , 0.08, "night", transform=ax[0].transAxes, verticalalignment='top', bbox=props)

        ax[0].set_ylabel("Height above ground level (m)",fontsize=14)
        ax[0].set_xlabel("θ (${^o}$C)",fontsize=14)


        ax[1].set_xticks(x_ticks)
        ax[1].set_xticklabels(x_ticklabels)
        ax[1].set_xlim(20,70)
        ax[1].set_ylim(0,2500)
        ax[1].text(0.02, 0.95, "(b)", transform=ax[1].transAxes, verticalalignment='top', bbox=props)
        ax[1].text(0.02, 0.08, "afternoon", transform=ax[1].transAxes, verticalalignment='top', bbox=props)
        ax[1].text(0.9 , 0.08, "night", transform=ax[1].transAxes, verticalalignment='top', bbox=props)
        ax[1].set_xlabel("θ (${^o}$C)",fontsize=14)
        # ax[1].set_ylabel("Geopotential Height (m)",fontsize=14)
        # end loc_num

        # =============== add legend ===============
        custom_lines = [Line2D([0], [0], color="orange", lw=1),
                        Line2D([0], [0], color="green", lw=1),
                        Line2D([0], [0], color="red", lw=1),
                        Line2D([0], [0], color="blue", lw=1)]

        fig.legend(custom_lines, ['pre-hw FD', 'pre-hw GW', 'hw FD', 'hw GW'],
                   loc='upper right', bbox_to_anchor=(0.83, 0.95), frameon=False)
        # fig.tight_layout()

        # ============ savefig ============
        message = 'theta_profile_' +case_name

        plt.savefig('./plots/figures/'+message+'.png',dpi=300)

if __name__ == "__main__":

    # #######################
    #        Set decks      #
    # #######################

    path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"

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


    plot_theta_profile(path, case_names, periods, time_ss, time_es,
                loc_lats=loc_lats, loc_lons=loc_lons)
