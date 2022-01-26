#!/usr/bin/python

'''
Plot spitial map of land diagnosis and parameters from LIS-CABLE
1. per time step
2. time period average
'''

from netCDF4 import Dataset
import numpy as np
import pandas as pd
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from pygam import LinearGAM
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_range_diff
from common_utils import *

def mask_by_hw_thrd(time_s,Var):

    print("In mask_by_hw_thrd Var")

    AWAP_mask  = "/g/data/w35/mm3972/scripts/Heatwave/coupled_GW_HW/nc_file/AWAP_tmean_above_Tmn90th_1970-2019_WRF_domain.nc"
    awap_file  = Dataset(AWAP_mask, mode='r')
    tmean      = awap_file.variables["tmean"][:]

    AWAP_time  = nc.num2date(awap_file.variables['time'][:],awap_file.variables['time'].units,
                 only_use_cftime_datetimes=False, only_use_python_datetimes=True)

    var        = Var

    for i in np.arange(len(Var[:,0,0])):
        date   = time_s + timedelta(days=int(i))

        # print("date ", date)
        # print("AWAP_time[AWAP_time == date]", AWAP_time[AWAP_time == date])

        var[i,:,:] = np.where(np.isnan(tmean[AWAP_time == date, :,:]),  np.nan, Var[i,:,:]) # awap_t == current for awap_t in AWAP_time

    # print(var)

    return var

def get_x_y_clr_val(df,var_name,cnt):

    iveg_num   = [2, 9, 5, 6, 14]

    normal_var = ["tair","tmax","tmin","qle","qh","tair","vegt","fw","sm"]
    et_var     = ["tveg","evap"]

    if var_name == "wtd":
        var = df[var_name][df.iveg == iveg_num[cnt]].values/1000.
    elif var_name in et_var:
        var = df[var_name][df.iveg == iveg_num[cnt]].values*3600*24
    else:
        var = df[var_name][df.iveg == iveg_num[cnt]].values

    return var

def plot_spatial_land(file_paths, time_s, time_e, seconds=None,loc_lat=None, loc_lon=None, message=None):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    file1  = Dataset(file_paths[0], mode='r')

    Time   = nc.num2date(file1.variables['time'][:],file1.variables['time'].units,
            only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time   = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)

    Tair1  = file1.variables['Tair_f_inst']
    tair1  = spital_var(time,Tair1,time_s,time_e,seconds)
    tmax1  = spital_var_max(time,Tair1,time_s,time_e,seconds)
    tmin1  = spital_var_min(time,Tair1,time_s,time_e,seconds)

    Qle1   = file1.variables['Qle_tavg']
    qle1   = spital_var(time,Qle1,time_s,time_e,seconds)
    Qh1    = file1.variables['Qh_tavg']
    qh1    = spital_var(time,Qh1,time_s,time_e,seconds)
    TVeg1  = file1.variables['TVeg_tavg']
    tveg1  = spital_var(time,TVeg1,time_s,time_e,seconds)
    Evap1  = file1.variables['Evap_tavg']
    evap1  = spital_var(time,Evap1,time_s,time_e,seconds)
    VegT1  = file1.variables['VegT_tavg']
    vegt1  = spital_var(time,VegT1,time_s,time_e,seconds)
    FW1    = file1.variables['FWsoil_tavg']
    fw1    = spital_var(time,FW1,time_s,time_e,seconds)

    if len(file_paths) > 1:
        file2  = Dataset(file_paths[1], mode='r')
        Tair2  = file2.variables['Tair_f_inst']

        tair2  = spital_var(time,Tair2,time_s,time_e,seconds)
        tmax2  = spital_var_max(time,Tair2,time_s,time_e,seconds)
        tmin2  = spital_var_min(time,Tair2,time_s,time_e,seconds)

        Qle2   = file2.variables['Qle_tavg']
        qle2   = spital_var(time,Qle2,time_s,time_e,seconds)
        Qh2    = file2.variables['Qh_tavg']
        qh2    = spital_var(time,Qh2,time_s,time_e,seconds)
        TVeg2  = file2.variables['TVeg_tavg']
        tveg2  = spital_var(time,TVeg2,time_s,time_e,seconds)
        Evap2  = file2.variables['Evap_tavg']
        evap2  = spital_var(time,Evap2,time_s,time_e,seconds)
        VegT2  = file2.variables['VegT_tavg']
        vegt2  = spital_var(time,VegT2,time_s,time_e,seconds)
        FW2    = file2.variables['FWsoil_tavg']
        fw2    = spital_var(time,FW2,time_s,time_e,seconds)

        iveg    = file2.variables["Landcover_inst"][0,:,:]
        WTD_tmp = file2.variables["WaterTableD_tavg"]
        wtd     = spital_var(time,WTD_tmp,time_s,time_e,seconds)

        tair     = tair2 - tair1
        tmax     = tmax2 - tmax1
        tmin     = tmin2 - tmin1
        qle      = qle2 - qle1
        qh       = qh2 - qh1
        tveg     = tveg2 - tveg1
        evap     = evap2 - evap1
        vegt     = vegt2 - vegt1
        fw       = fw2 - fw1
    else:
        tair     = tair1
        tmax     = tmax1
        tmin     = tmin1
        qle      = qle1
        qh       = qh1
        tveg     = tveg1
        evap     = evap1
        vegt     = vegt1
        fw       = fw1


    # # remove ocean and boundary pixels
    # nlon = len(tair[0,6:-5])
    # nlat = len(tair[6:-5,0])

    # the new array doesn't set a missing value, and the missing value from old array will become nan
    # iveg_shrink    = np.zeros(nlat * nlon)

    iveg_shrink    = np.reshape(iveg[6:-5,6:-5],-1)
    wtd_shrink     = np.reshape(wtd[6:-5,6:-5],-1)
    tair_shrink    = np.reshape(tair[6:-5,6:-5],-1)
    tmax_shrink    = np.reshape(tmax[6:-5,6:-5],-1)
    tmin_shrink    = np.reshape(tmin[6:-5,6:-5],-1)
    qle_shrink     = np.reshape(qle[6:-5,6:-5],-1)
    qh_shrink      = np.reshape(qh[6:-5,6:-5],-1)
    tveg_shrink    = np.reshape(tveg[6:-5,6:-5],-1)
    evap_shrink    = np.reshape(evap[6:-5,6:-5],-1)
    vegt_shrink    = np.reshape(vegt[6:-5,6:-5],-1)
    fw_shrink      = np.reshape(fw[6:-5,6:-5],-1)

    # print("wtd_shrink",~ np.isnan(wtd_shrink))
    iveg_1D    = iveg_shrink[~ np.isnan(wtd_shrink)]
    wtd_1D     = wtd_shrink[~ np.isnan(wtd_shrink)]
    tair_1D    = tair_shrink[~ np.isnan(wtd_shrink)]
    tmax_1D    = tmax_shrink[~ np.isnan(wtd_shrink)]
    tmin_1D    = tmin_shrink[~ np.isnan(wtd_shrink)]
    qle_1D     = qle_shrink[~ np.isnan(wtd_shrink)]
    qh_1D      = qh_shrink[~ np.isnan(wtd_shrink)]
    tveg_1D    = tveg_shrink[~ np.isnan(wtd_shrink)]
    evap_1D    = evap_shrink[~ np.isnan(wtd_shrink)]
    vegt_1D    = vegt_shrink[~ np.isnan(wtd_shrink)]
    fw_1D      = fw_shrink[~ np.isnan(wtd_shrink)]

    df         = pd.DataFrame(iveg_1D, columns=['iveg'])
    df['wtd']  = wtd_1D
    df['tair'] = tair_1D
    df['tmax'] = tmax_1D
    df['tmin'] = tmin_1D
    df['qle']  = qle_1D
    df['qh']   = qh_1D
    df['tveg'] = tveg_1D
    df['evap'] = evap_1D
    df['vegt'] = vegt_1D
    df['fw']   = fw_1D

    # df['rank']     = df[var_name].rank()
    # df_sort        = df.sort_values(by=var_name)

    # ============ Preset of plotting ============
    # cmap     = cm.Set2(np.arange(5))
    pft      = ["BEF","shrub","grass","crop","barren"]
    order    = ["(a)","(b)","(c)","(d)","(e)"]
    iveg_num = [2, 5, 6, 9, 14]
    cmap     = cm.Set2(np.arange(5))
    markers  = ["o","o","o","o","o"]
    mrk_sz   = 4


    # ============ Plotting ============
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[12,8], sharex=True, sharey=True, squeeze=True)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
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

    cnt = 0

    for i in np.arange(2):
        for j in np.arange(2):
            if cnt < 5:
                wtd_val    = df['wtd'][df.iveg == iveg_num[cnt]].values/1000.
                tair_val   = df['tair'][df.iveg == iveg_num[cnt]].values
                tmax_val   = df['tmax'][df.iveg == iveg_num[cnt]].values
                tmin_val   = df['tmin'][df.iveg == iveg_num[cnt]].values
                qle_val    = df['qle'][df.iveg == iveg_num[cnt]].values
                qh_val     = df['qh'][df.iveg == iveg_num[cnt]].values
                tveg_val   = df['tveg'][df.iveg == iveg_num[cnt]].values*3600*24
                evap_val   = df['evap'][df.iveg == iveg_num[cnt]].values*3600*24
                vegt_val   = df['vegt'][df.iveg == iveg_num[cnt]].values
                fw_val     = df['fw'][df.iveg == iveg_num[cnt]].values

                xxx_value  = wtd_val #fw_val #wtd_val
                var_value  = tmax_val
                clr_value  = fw_val #wtd_val #fw_val
                val_min    = -0.1 #0.
                val_max    = 0.8  #15.

                ax[i,j].set_xlim(-0.5,18.5)
                # ax[i,j].set_xlim(-0.05,0.8)
                ax[i,j].set_ylim(-2.49,0.6)

                if i == 1:
                    ax[i,j].set_xlabel("WTD (m)", fontsize=12)
                    # ax[i,j].set_xlabel("Δ$β$ (-)", fontsize=12)
                if j == 0:
                    ax[i,j].set_ylabel("ΔT$_{max}$ ($\mathregular{^o}$C)", fontsize=12)

                # gam fitting
                # gam    = LinearGAM(n_splines=4).gridsearch(xxx_value[:,None], var_value) # n_splines=22
                # x_pred = np.linspace(min(xxx_value), max(xxx_value), num=100)
                # y_pred = gam.predict(x_pred)
                # y_int  = gam.confidence_intervals(x_pred, width=.95)

                sct    = ax[i,j].scatter(xxx_value, var_value, s=mrk_sz, vmin=val_min, vmax=val_max, marker=markers[cnt],c=clr_value, cmap="BrBG", label=pft[cnt]) # alpha=0.7,
                ax[i,j].axhline(y=0, color="gray", linestyle='--')
                if cnt == 3:
                    cbar   = plt.colorbar(sct, ax=ax, ticklocation="right", pad=0.01, orientation="vertical", aspect=40, shrink=1.)
                        #  norm=colors.NoNorm(),boundaries=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7], ticks=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7]) # cax=cax,
                    cbar.set_label('Δ$β$ (-)',  loc='center') # rotation=270,
                    # cbar.set_label('WTD (m)',  loc='center') # rotation=270,
                ax[i,j].text(0.02, 0.95, order[cnt]+" "+pft[cnt], transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)


                # plt.title(pft[i], size=12)
                # clb.set_clim(0, 1)

                # ax[i,j].plot(x_pred, y_pred, color="black", ls='-', lw=2.0, zorder=10)
                # ax[i,j].fill_between(x_pred, y_int[:, 0], y_int[:, 1], alpha=0.2,
                #                 facecolor='grey', zorder=10)

                # ax[i,j].legend()

                cnt = cnt + 1

    fig.savefig("./plots/figures/scatter_delta_Tmax_WTD_PFT_"+message,bbox_inches='tight')#, pad_inches=0.1)

def plot_spatial_land_days(case_names, periods, time_ss, time_es, message=None):

    hw_mask = True

    # ============ Setting for plotting ============
    # # wtd-tmax-fw
    # x_y_clr  = ["wtd","tmax","fw"]
    # labels   = ["WTD (m)",           # x-axis
    #             "ΔT$_{max}$ ($\mathregular{^o}$C)", # y-axis
    #             "Δ$β$ (-)"]          # colorbar
    # x_min   = 0
    # x_max   = 21
    # y_min   = -6
    # y_max   = 3.
    # clr_min = -0.1#-0.8
    # clr_max = 0.8

    # # wtd-tmax-tmax_gw
    # x_y_clr  = ["wtd","tmax","tmax_gw"]
    # labels   = ["WTD (m)",          # x-axis
    #             "ΔT$_{max}$ ($\mathregular{^o}$C)", # y-axis
    #             "T$_{gw,max}$ ($\mathregular{^o}$C)"]          # colorbar
    # x_min  = 0
    # x_max  = 21
    # y_min  = -6.
    # y_max  = 3.
    # clr_min = 30 #-0.8
    # clr_max = 50
    
    # tmax_gw-tmax-fw
    x_y_clr  = ["tmax_gw","tmax","fw"]
    labels   = ["T$_{gw,max}$ ($\mathregular{^o}$C)",          # x-axis
                "ΔT$_{max}$ ($\mathregular{^o}$C)", # y-axis
                "Δ$β$ (-)"]          # colorbar
    x_min  = 15
    x_max  = 52
    y_min  = -6
    y_max  = 3
    clr_min = -0.1 #-0.8
    clr_max = 0.8
    
    # wtd-sm-fw
    # x_y_clr  = ["wtd","sm","fw"]
    # labels   = ["WTD (m)",           # x-axis
    #             "ΔSM (m${^3}$ m$^{-3}$)", # y-axis
    #             "Δ$β$ (-)"]          # colorbar
    # x_min  = -0.5
    # x_max  = 20
    # y_min  = -0.01
    # y_max  = 0.2
    # clr_min = -0.1#-0.8
    # clr_max = 0.8

    # # fw-sm-wtd
    # x_y_clr  = ["fw","sm","wtd"]
    # labels   = ["Δ$β$ (-)",           # x-axis
    #             "ΔSM (m${^3}$ m$^{-3}$)", # y-axis
    #             "WTD (m)"]          # colorbar
    # x_min  = -0.5
    # x_max  = 20
    # y_min  = -0.01
    # y_max  = 0.2
    # clr_min = -0.1#-0.8
    # clr_max = 0.8

    # sm-fw-wtd
    # x_y_clr  = ["sm","fw","wtd"]
    # labels   = ["ΔSM (m${^3}$ m$^{-3}$)","Δ$β$ (-)",
    #             "WTD (m)"]          # colorbar
    # x_min  = -0.01
    # x_max  = 0.2
    # y_min  = -0.2
    # y_max  = 1.
    # clr_min = -0.5#-0.8
    # clr_max = 20.

    # # fw-tmax-tmax_gw
    # x_y_clr  = ["fw","tmax","tmax_gw"]
    # labels   = ["Δ$β$ (-)",           # x-axis
    #             "ΔT$_{max}$ ($\mathregular{^o}$C)", # y-axis
    #             "T$_{gw,max}$ ($\mathregular{^o}$C)"]          # colorbar
    # x_min  = -0.2
    # x_max  = 1.0
    # y_min  = -6.
    # y_max  = 3.
    # clr_min = 15 #-0.8
    # clr_max = 50

    # ============ Setting for plotting ============
    pft      = ["BEF","crop","shrub","grass","barren"]
    order    = ["(a)","(b)","(c)","(d)","(e)"]
    iveg_num = [2, 9, 5, 6, 14]
    cmap     = plt.cm.BrBG #YlOrBr #coolwarm_r
    cmap_new = truncate_colormap(cmap, minval=0.5-0.5/8, maxval=1.)
    markers  = ["o","o","o","^","s"]
    mrk_sz   = 1.5

    fig, ax  = plt.subplots(nrows=2, ncols=2, figsize=[12,8],sharex=True, sharey=True, squeeze=True) #
    plt.subplots_adjust(wspace=0.0, hspace=0.0)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
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

    # ============ Read data ============
    for i,case_name in enumerate(case_names):

        time_s = time_ss[i]
        time_e = time_es[i]
        period = periods[i]

        cpl_land_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'

        cpl_land_file_gw  = cpl_land_file + '/LIS.CABLE.'+period+'_gw.nc'  # land output of wrf-cable run
        cpl_land_file_fd  = cpl_land_file + '/LIS.CABLE.'+period+'_fd.nc'  # land output of wrf-cable run

        file_paths        = [cpl_land_file_fd,cpl_land_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw

        # Open the NetCDF4 file (add a directory path if necessary) for reading:
        file1  = Dataset(file_paths[0], mode='r')

        Time   = nc.num2date(file1.variables['time'][:],file1.variables['time'].units,
                 only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        time   = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)

        Tair1  = file1.variables['Tair_f_inst']
        tair1  = time_clip_to_day(time,Tair1,time_s,time_e)
        tmax1  = time_clip_to_day_max(time,Tair1,time_s,time_e)
        tmin1  = time_clip_to_day_min(time,Tair1,time_s,time_e)

        Qle1   = file1.variables['Qle_tavg']
        qle1   = time_clip_to_day(time,Qle1,time_s,time_e)
        Qh1    = file1.variables['Qh_tavg']
        qh1    = time_clip_to_day(time,Qh1,time_s,time_e)
        TVeg1  = file1.variables['TVeg_tavg']
        tveg1  = time_clip_to_day(time,TVeg1,time_s,time_e)
        Evap1  = file1.variables['Evap_tavg']
        evap1  = time_clip_to_day(time,Evap1,time_s,time_e)
        VegT1  = file1.variables['VegT_tavg']
        vegt1  = time_clip_to_day(time,VegT1,time_s,time_e)
        FW1    = file1.variables['FWsoil_tavg']
        fw1    = time_clip_to_day(time,FW1,time_s,time_e)

        # top 64cm SM
        SM1    = ((file1.variables['SoilMoist_inst'][:,0,:,:]*0.022 +
                   file1.variables['SoilMoist_inst'][:,1,:,:]*0.058 +
                   file1.variables['SoilMoist_inst'][:,2,:,:]*0.154 +
                   file1.variables['SoilMoist_inst'][:,3,:,:]*0.409)/
                  (0.022+0.058+0.154+0.409))
        sm1    = time_clip_to_day(time,SM1,time_s,time_e)

        if len(file_paths) > 1:
            file2  = Dataset(file_paths[1], mode='r')
            Tair2  = file2.variables['Tair_f_inst']

            tair2  = time_clip_to_day(time,Tair2,time_s,time_e)
            tmax2  = time_clip_to_day_max(time,Tair2,time_s,time_e)
            tmin2  = time_clip_to_day_min(time,Tair2,time_s,time_e)

            Qle2   = file2.variables['Qle_tavg']
            qle2   = time_clip_to_day(time,Qle2,time_s,time_e)
            Qh2    = file2.variables['Qh_tavg']
            qh2    = time_clip_to_day(time,Qh2,time_s,time_e)
            TVeg2  = file2.variables['TVeg_tavg']
            tveg2  = time_clip_to_day(time,TVeg2,time_s,time_e)
            Evap2  = file2.variables['Evap_tavg']
            evap2  = time_clip_to_day(time,Evap2,time_s,time_e)
            VegT2  = file2.variables['VegT_tavg']
            vegt2  = time_clip_to_day(time,VegT2,time_s,time_e)
            FW2    = file2.variables['FWsoil_tavg']
            fw2    = time_clip_to_day(time,FW2,time_s,time_e)

            # top 64cm SM
            SM2    = ((file2.variables['SoilMoist_inst'][:,0,:,:]*0.022 +
                    file2.variables['SoilMoist_inst'][:,1,:,:]*0.058 +
                    file2.variables['SoilMoist_inst'][:,2,:,:]*0.154 +
                    file2.variables['SoilMoist_inst'][:,3,:,:]*0.409)/
                    (0.022+0.058+0.154+0.409))
            sm2    = time_clip_to_day(time,SM2,time_s,time_e)

            Iveg    = file2.variables["Landcover_inst"]
            iveg    = time_clip_to_day(time,Iveg,time_s,time_e)
            WTD     = file2.variables["WaterTableD_tavg"]
            wtd     = time_clip_to_day(time,WTD,time_s,time_e)

            tair     = tair2 - tair1
            tmax     = tmax2 - tmax1
            tmin     = tmin2 - tmin1
            qle      = qle2 - qle1
            qh       = qh2 - qh1
            tveg     = tveg2 - tveg1
            evap     = evap2 - evap1
            vegt     = vegt2 - vegt1
            fw       = fw2 - fw1
            sm       = sm2 - sm1 ###
        else:
            tair     = tair1
            tmax     = tmax1
            tmin     = tmin1
            qle      = qle1
            qh       = qh1
            tveg     = tveg1
            evap     = evap1
            vegt     = vegt1
            fw       = fw1

        if hw_mask:
            # use wtd to mask other arrays
            # print("np.shape(wtd)", np.shape(wtd))
            wtd     = mask_by_hw_thrd(time_s,wtd)
            print("np.shape(wtd)", np.shape(wtd))

        iveg_shrink    = np.reshape(iveg[:,6:-5,6:-5],-1)
        wtd_shrink     = np.reshape(wtd[:,6:-5,6:-5],-1)
        tair_shrink    = np.reshape(tair[:,6:-5,6:-5],-1)
        tmax_shrink    = np.reshape(tmax[:,6:-5,6:-5],-1)
        tmax_gw_shrink = np.reshape(tmax2[:,6:-5,6:-5],-1)
        tmin_shrink    = np.reshape(tmin[:,6:-5,6:-5],-1)
        qle_shrink     = np.reshape(qle[:,6:-5,6:-5],-1)
        qh_shrink      = np.reshape(qh[:,6:-5,6:-5],-1)
        tveg_shrink    = np.reshape(tveg[:,6:-5,6:-5],-1)
        evap_shrink    = np.reshape(evap[:,6:-5,6:-5],-1)
        vegt_shrink    = np.reshape(vegt[:,6:-5,6:-5],-1)
        fw_shrink      = np.reshape(fw[:,6:-5,6:-5],-1)
        sm_shrink      = np.reshape(sm[:,6:-5,6:-5],-1)

        # print(np.shape(iveg_shrink),np.shape(wtd_shrink),np.shape(tair_shrink),np.shape(tmax_shrink),np.shape(tmax_gw_shrink),
        #       np.shape(tmin_shrink),np.shape(qle_shrink),np.shape(qh_shrink),np.shape(tveg_shrink),np.shape(evap_shrink),
        #       np.shape(vegt_shrink),np.shape(fw_shrink),np.shape(sm_shrink))

        # mask out ocean pixels
        iveg_1D    = iveg_shrink[~ np.isnan(wtd_shrink)]
        wtd_1D     = wtd_shrink[~ np.isnan(wtd_shrink)]
        tair_1D    = tair_shrink[~ np.isnan(wtd_shrink)]
        tmax_1D    = tmax_shrink[~ np.isnan(wtd_shrink)]
        tmax_gw_1D = tmax_gw_shrink[~ np.isnan(wtd_shrink)]
        tmin_1D    = tmin_shrink[~ np.isnan(wtd_shrink)]
        qle_1D     = qle_shrink[~ np.isnan(wtd_shrink)]
        qh_1D      = qh_shrink[~ np.isnan(wtd_shrink)]
        tveg_1D    = tveg_shrink[~ np.isnan(wtd_shrink)]
        evap_1D    = evap_shrink[~ np.isnan(wtd_shrink)]
        vegt_1D    = vegt_shrink[~ np.isnan(wtd_shrink)]
        fw_1D      = fw_shrink[~ np.isnan(wtd_shrink)]
        sm_1D      = sm_shrink[~ np.isnan(wtd_shrink)]

        # print(sm_1D)

        if i == 0:
            df         = pd.DataFrame(iveg_1D, columns=['iveg'])
            df['hw']   = case_name
            df['wtd']  = wtd_1D
            df['tair'] = tair_1D
            df['tmax'] = tmax_1D
            df['tmax_gw'] = tmax_gw_1D - 273.15
            df['tmin'] = tmin_1D
            df['qle']  = qle_1D
            df['qh']   = qh_1D
            df['tveg'] = tveg_1D
            df['evap'] = evap_1D
            df['vegt'] = vegt_1D
            df['fw']   = fw_1D
            df['sm']   = sm_1D
            print(df)
        else:
            df2         = pd.DataFrame(iveg_1D, columns=['iveg'])
            df2['hw']   = case_name
            df2['wtd']  = wtd_1D
            df2['tair'] = tair_1D
            df2['tmax'] = tmax_1D
            df2['tmax_gw'] = tmax_gw_1D - 273.15
            df2['tmin'] = tmin_1D
            df2['qle']  = qle_1D
            df2['qh']   = qh_1D
            df2['tveg'] = tveg_1D
            df2['evap'] = evap_1D
            df2['vegt'] = vegt_1D
            df2['fw']   = fw_1D
            df2['sm']   = sm_1D
            df = df.append(df2)
            print(df2)
            df2 = None
    print(df)

    # ============ Setting for plotting ============
    cnt = 0

    for j in np.arange(2):
        # if j == 1:
        #     y_min  = -2.2
        #     y_max  = 2.2
        for k in np.arange(2):

            xxx_value  = get_x_y_clr_val(df,x_y_clr[0],cnt)
            var_value  = get_x_y_clr_val(df,x_y_clr[1],cnt)
            clr_value  = get_x_y_clr_val(df,x_y_clr[2],cnt)
            median_val = np.nanmedian(xxx_value)
            median_var_val = np.nanmedian(var_value)

            ax[j,k].set_xlim(x_min,x_max)
            ax[j,k].set_ylim(y_min,y_max)

            if j == 1:
                # i == 0 and
                ax[j,k].set_xlabel(labels[0], fontsize=12)
                # ax[j,k].set_xlabel("Δ$β$ (-)", fontsize=12)
            if k == 0:
                # i == 0 and
                ax[j,k].set_ylabel(labels[1], fontsize=12)

            # gam fitting
            # gam    = LinearGAM(n_splines=4).gridsearch(xxx_value[:,None], var_value) # n_splines=22
            # x_pred = np.linspace(min(xxx_value), max(xxx_value), num=100)
            # y_pred = gam.predict(x_pred)
            # y_int  = gam.confidence_intervals(x_pred, width=.95)

            sct    = ax[j,k].scatter(xxx_value, var_value, s=mrk_sz, vmin=clr_min, vmax=clr_max, marker=markers[0],
                                    c=clr_value, cmap=cmap, label=pft[cnt]) # ,facecolors='none' marker=markers[i]  alpha=0., edgecolors=clr_value,

            # if i == 0:
            ax[j,k].axvline(x=median_val, color="gray", linestyle='--')
            ax[j,k].axhline(y=median_var_val, color="gray", linestyle='--')
            ax[j,k].text(0.02, 0.95, order[cnt]+" "+pft[cnt], transform=ax[j,k].transAxes, fontsize=14, verticalalignment='top', bbox=props)

            ax[j,k].text(0.4, 0.18, "Median", transform=ax[j,k].transAxes, fontsize=12, verticalalignment='top', bbox=props)
            ax[j,k].text(0.4, 0.12, "WTD     "+str("{:.2f}".format(median_val))+" m", transform=ax[j,k].transAxes, fontsize=12, verticalalignment='top', bbox=props)
            ax[j,k].text(0.4, 0.07, "ΔT$_{max}$  "+str("{:.2f}".format(median_var_val))+" $\mathregular{^o}$C",
                         transform=ax[j,k].transAxes, fontsize=12, verticalalignment='top', bbox=props)

            if cnt == 3:
                cbar   = plt.colorbar(sct, ax=ax, ticklocation="right", pad=0.01, orientation="vertical",
                                        aspect=40, shrink=1.,extend="both") # drawedges=True,
                    #  norm=colors.NoNorm(),boundaries=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7], ticks=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7]) # cax=cax,
                cbar.set_label(labels[2],  loc='center') # rotation=270

            cnt = cnt + 1

            # plt.title(pft[j], size=12)
            # clb.set_clim(0, 1)

            # ax[j,k].plot(x_pred, y_pred, color="black", ls='-', lw=2.0, zorder=10)
            # ax[j,k].fill_between(x_pred, y_int[:, 0], y_int[:, 1], alpha=0.2,
            #                 facecolor='grey', zorder=10)

        # df = None
    message = x_y_clr[0]+"_"+x_y_clr[1]+"_"+x_y_clr[2]
    fig.savefig("./plots/figures/scatter_"+message+"_days",bbox_inches='tight')#, pad_inches=0.1)

def plot_spatial_tair_fwsoil(file_paths, time_s, time_e, seconds=None,loc_lat=None, loc_lon=None, message=None):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    file1  = Dataset(file_paths[0], mode='r')

    Time   = nc.num2date(file1.variables['time'][:],file1.variables['time'].units,
            only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time   = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)

    Tair1  = file1.variables['Tair_f_inst']
    tair1  = spital_var(time,Tair1,time_s,time_e,seconds)
    tmax1  = spital_var_max(time,Tair1,time_s,time_e,seconds)
    tmin1  = spital_var_min(time,Tair1,time_s,time_e,seconds)

    Qle1   = file1.variables['Qle_tavg']
    qle1   = spital_var(time,Qle1,time_s,time_e,seconds)
    Qh1    = file1.variables['Qh_tavg']
    qh1    = spital_var(time,Qh1,time_s,time_e,seconds)
    TVeg1  = file1.variables['TVeg_tavg']
    tveg1  = spital_var(time,TVeg1,time_s,time_e,seconds)
    Evap1  = file1.variables['Evap_tavg']
    evap1  = spital_var(time,Evap1,time_s,time_e,seconds)
    VegT1  = file1.variables['VegT_tavg']
    vegt1  = spital_var(time,VegT1,time_s,time_e,seconds)
    FW1    = file1.variables['FWsoil_tavg']
    fw1    = spital_var(time,FW1,time_s,time_e,seconds)

    if len(file_paths) > 1:
        file2  = Dataset(file_paths[1], mode='r')
        Tair2  = file2.variables['Tair_f_inst']

        tair2  = spital_var(time,Tair2,time_s,time_e,seconds)
        tmax2  = spital_var_max(time,Tair2,time_s,time_e,seconds)
        tmin2  = spital_var_min(time,Tair2,time_s,time_e,seconds)

        Qle2   = file2.variables['Qle_tavg']
        qle2   = spital_var(time,Qle2,time_s,time_e,seconds)
        Qh2    = file2.variables['Qh_tavg']
        qh2    = spital_var(time,Qh2,time_s,time_e,seconds)
        TVeg2  = file2.variables['TVeg_tavg']
        tveg2  = spital_var(time,TVeg2,time_s,time_e,seconds)
        Evap2  = file2.variables['Evap_tavg']
        evap2  = spital_var(time,Evap2,time_s,time_e,seconds)
        VegT2  = file2.variables['VegT_tavg']
        vegt2  = spital_var(time,VegT2,time_s,time_e,seconds)
        FW2    = file2.variables['FWsoil_tavg']
        fw2    = spital_var(time,FW2,time_s,time_e,seconds)

        iveg    = file2.variables["Landcover_inst"][0,:,:]
        WTD_tmp = file2.variables["WaterTableD_tavg"]
        wtd     = spital_var(time,WTD_tmp,time_s,time_e,seconds)

        tair     = tair2  -273.15 #- tair1
        tmax     = tmax2 - tmax1
        tmin     = tmin2 - tmin1
        qle      = qle2 - qle1
        qh       = qh2 - qh1
        tveg     = tveg2 - tveg1
        evap     = evap2 - evap1
        vegt     = vegt2 - vegt1
        fw       = fw2 - fw1
    else:
        tair     = tair1
        tmax     = tmax1
        tmin     = tmin1
        qle      = qle1
        qh       = qh1
        tveg     = tveg1
        evap     = evap1
        vegt     = vegt1
        fw       = fw1


    # # remove ocean and boundary pixels
    # nlon = len(tair[0,6:-5])
    # nlat = len(tair[6:-5,0])

    # the new array doesn't set a missing value, and the missing value from old array will become nan
    # iveg_shrink    = np.zeros(nlat * nlon)

    iveg_shrink    = np.reshape(iveg[6:-5,6:-5],-1)
    wtd_shrink     = np.reshape(wtd[6:-5,6:-5],-1)
    tair_shrink    = np.reshape(tair[6:-5,6:-5],-1)
    tmax_shrink    = np.reshape(tmax[6:-5,6:-5],-1)
    tmin_shrink    = np.reshape(tmin[6:-5,6:-5],-1)
    qle_shrink     = np.reshape(qle[6:-5,6:-5],-1)
    qh_shrink      = np.reshape(qh[6:-5,6:-5],-1)
    tveg_shrink    = np.reshape(tveg[6:-5,6:-5],-1)
    evap_shrink    = np.reshape(evap[6:-5,6:-5],-1)
    vegt_shrink    = np.reshape(vegt[6:-5,6:-5],-1)
    fw_shrink      = np.reshape(fw[6:-5,6:-5],-1)

    # print("wtd_shrink",~ np.isnan(wtd_shrink))
    iveg_1D    = iveg_shrink[~ np.isnan(wtd_shrink)]
    wtd_1D     = wtd_shrink[~ np.isnan(wtd_shrink)]
    tair_1D    = tair_shrink[~ np.isnan(wtd_shrink)]
    tmax_1D    = tmax_shrink[~ np.isnan(wtd_shrink)]
    tmin_1D    = tmin_shrink[~ np.isnan(wtd_shrink)]
    qle_1D     = qle_shrink[~ np.isnan(wtd_shrink)]
    qh_1D      = qh_shrink[~ np.isnan(wtd_shrink)]
    tveg_1D    = tveg_shrink[~ np.isnan(wtd_shrink)]
    evap_1D    = evap_shrink[~ np.isnan(wtd_shrink)]
    vegt_1D    = vegt_shrink[~ np.isnan(wtd_shrink)]
    fw_1D      = fw_shrink[~ np.isnan(wtd_shrink)]

    df         = pd.DataFrame(iveg_1D, columns=['iveg'])
    df['wtd']  = wtd_1D
    df['tair'] = tair_1D
    df['tmax'] = tmax_1D
    df['tmin'] = tmin_1D
    df['qle']  = qle_1D
    df['qh']   = qh_1D
    df['tveg'] = tveg_1D
    df['evap'] = evap_1D
    df['vegt'] = vegt_1D
    df['fw']   = fw_1D

    # df['rank']     = df[var_name].rank()
    # df_sort        = df.sort_values(by=var_name)

    # ============ Preset of plotting ============
    # cmap     = cm.Set2(np.arange(5))
    pft      = ["BEF","shrub","grass","crop","barren"]
    order    = ["(a)","(b)","(c)","(d)","(e)"]
    iveg_num = [2, 5, 6, 9, 14]
    cmap     = cm.Set2(np.arange(5))
    markers  = ["o","o","o","o","o"]
    mrk_sz   = 4


    # ============ Plotting ============
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[12,8],sharex=True, sharey=True, squeeze=True)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
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

    cnt = 0

    for i in np.arange(2):
        for j in np.arange(2):
            if cnt < 5:
                wtd_val    = df['wtd'][df.iveg == iveg_num[cnt]].values/1000.
                tair_val   = df['tair'][df.iveg == iveg_num[cnt]].values
                tmax_val   = df['tmax'][df.iveg == iveg_num[cnt]].values
                tmin_val   = df['tmin'][df.iveg == iveg_num[cnt]].values
                qle_val    = df['qle'][df.iveg == iveg_num[cnt]].values
                qh_val     = df['qh'][df.iveg == iveg_num[cnt]].values
                tveg_val   = df['tveg'][df.iveg == iveg_num[cnt]].values*3600*24
                evap_val   = df['evap'][df.iveg == iveg_num[cnt]].values*3600*24
                vegt_val   = df['vegt'][df.iveg == iveg_num[cnt]].values
                fw_val     = df['fw'][df.iveg == iveg_num[cnt]].values

                xxx_value  = tair_val #fw_val #wtd_val
                var_value  = tmax_val
                clr_value  = fw_val #wtd_val #fw_val
                val_min    = -0.1 #0.
                val_max    = 0.8  #15.

                ax[i,j].set_xlim(15,45)
                # ax[i,j].set_xlim(-0.05,0.8)
                ax[i,j].set_ylim(-2.49,0.6)

                if i == 1:
                    ax[i,j].set_xlabel("T ($\mathregular{^o}$C)", fontsize=12)
                    # ax[i,j].set_xlabel("Δ$β$ (-)", fontsize=12)
                if j == 0:
                    ax[i,j].set_ylabel("ΔTmax ($\mathregular{^o}$C)", fontsize=12)

                # gam fitting
                gam    = LinearGAM(n_splines=4).gridsearch(xxx_value[:,None], var_value) # n_splines=22
                x_pred = np.linspace(min(xxx_value), max(xxx_value), num=100)
                y_pred = gam.predict(x_pred)
                y_int  = gam.confidence_intervals(x_pred, width=.95)

                sct    = ax[i,j].scatter(xxx_value, var_value, s=mrk_sz, vmin=val_min, vmax=val_max, marker=markers[cnt],c=clr_value, cmap="BrBG", label=pft[cnt]) # alpha=0.7,
                ax[i,j].axhline(y=0, color="gray", linestyle='--')
                if cnt == 3:
                    cbar   = plt.colorbar(sct, ax=ax, ticklocation="right", pad=0.01, orientation="vertical", aspect=40, shrink=1.)
                        #  norm=colors.NoNorm(),boundaries=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7], ticks=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7]) # cax=cax,
                    cbar.set_label('Δ$β$ (-)',  loc='center') # rotation=270,
                    # cbar.set_label('WTD (m)',  loc='center') # rotation=270,
                ax[i,j].text(0.02, 0.95, order[cnt]+" "+pft[cnt], transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)


                # plt.title(pft[i], size=12)
                # clb.set_clim(0, 1)

                ax[i,j].plot(x_pred, y_pred, color="black", ls='-', lw=2.0, zorder=10)
                ax[i,j].fill_between(x_pred, y_int[:, 0], y_int[:, 1], alpha=0.2,
                                facecolor='grey', zorder=10)

                # ax[i,j].legend()

                cnt = cnt + 1

    fig.savefig("./plots/figures/scatter_delta_Tmax_Tair_FWsoil_PFT_"+message,bbox_inches='tight')#, pad_inches=0.1)

def plot_spatial_fwsoil_tair(file_paths, time_s, time_e, seconds=None,loc_lat=None, loc_lon=None, message=None):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    file1  = Dataset(file_paths[0], mode='r')

    Time   = nc.num2date(file1.variables['time'][:],file1.variables['time'].units,
            only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time   = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)

    Tair1  = file1.variables['Tair_f_inst']
    tair1  = spital_var(time,Tair1,time_s,time_e,seconds)
    tmax1  = spital_var_max(time,Tair1,time_s,time_e,seconds)
    tmin1  = spital_var_min(time,Tair1,time_s,time_e,seconds)

    Qle1   = file1.variables['Qle_tavg']
    qle1   = spital_var(time,Qle1,time_s,time_e,seconds)
    Qh1    = file1.variables['Qh_tavg']
    qh1    = spital_var(time,Qh1,time_s,time_e,seconds)
    TVeg1  = file1.variables['TVeg_tavg']
    tveg1  = spital_var(time,TVeg1,time_s,time_e,seconds)
    Evap1  = file1.variables['Evap_tavg']
    evap1  = spital_var(time,Evap1,time_s,time_e,seconds)
    VegT1  = file1.variables['VegT_tavg']
    vegt1  = spital_var(time,VegT1,time_s,time_e,seconds)
    FW1    = file1.variables['FWsoil_tavg']
    fw1    = spital_var(time,FW1,time_s,time_e,seconds)

    if len(file_paths) > 1:
        file2  = Dataset(file_paths[1], mode='r')
        Tair2  = file2.variables['Tair_f_inst']

        tair2  = spital_var(time,Tair2,time_s,time_e,seconds)
        tmax2  = spital_var_max(time,Tair2,time_s,time_e,seconds)
        tmin2  = spital_var_min(time,Tair2,time_s,time_e,seconds)

        Qle2   = file2.variables['Qle_tavg']
        qle2   = spital_var(time,Qle2,time_s,time_e,seconds)
        Qh2    = file2.variables['Qh_tavg']
        qh2    = spital_var(time,Qh2,time_s,time_e,seconds)
        TVeg2  = file2.variables['TVeg_tavg']
        tveg2  = spital_var(time,TVeg2,time_s,time_e,seconds)
        Evap2  = file2.variables['Evap_tavg']
        evap2  = spital_var(time,Evap2,time_s,time_e,seconds)
        VegT2  = file2.variables['VegT_tavg']
        vegt2  = spital_var(time,VegT2,time_s,time_e,seconds)
        FW2    = file2.variables['FWsoil_tavg']
        fw2    = spital_var(time,FW2,time_s,time_e,seconds)

        iveg    = file2.variables["Landcover_inst"][0,:,:]
        WTD_tmp = file2.variables["WaterTableD_tavg"]
        wtd     = spital_var(time,WTD_tmp,time_s,time_e,seconds)

        tair     = tair2 -273.15 #- tair1
        tmax     = tmax2 - tmax1
        tmin     = tmin2 - tmin1
        qle      = qle2 - qle1
        qh       = qh2 - qh1
        tveg     = tveg2 - tveg1
        evap     = evap2 - evap1
        vegt     = vegt2 - vegt1
        fw       = fw2 - fw1
    else:
        tair     = tair1
        tmax     = tmax1
        tmin     = tmin1
        qle      = qle1
        qh       = qh1
        tveg     = tveg1
        evap     = evap1
        vegt     = vegt1
        fw       = fw1


    # # remove ocean and boundary pixels
    # nlon = len(tair[0,6:-5])
    # nlat = len(tair[6:-5,0])

    # the new array doesn't set a missing value, and the missing value from old array will become nan
    # iveg_shrink    = np.zeros(nlat * nlon)

    iveg_shrink    = np.reshape(iveg[6:-5,6:-5],-1)
    wtd_shrink     = np.reshape(wtd[6:-5,6:-5],-1)
    tair_shrink    = np.reshape(tair[6:-5,6:-5],-1)
    tmax_shrink    = np.reshape(tmax[6:-5,6:-5],-1)
    tmin_shrink    = np.reshape(tmin[6:-5,6:-5],-1)
    qle_shrink     = np.reshape(qle[6:-5,6:-5],-1)
    qh_shrink      = np.reshape(qh[6:-5,6:-5],-1)
    tveg_shrink    = np.reshape(tveg[6:-5,6:-5],-1)
    evap_shrink    = np.reshape(evap[6:-5,6:-5],-1)
    vegt_shrink    = np.reshape(vegt[6:-5,6:-5],-1)
    fw_shrink      = np.reshape(fw[6:-5,6:-5],-1)

    # print("wtd_shrink",~ np.isnan(wtd_shrink))
    iveg_1D    = iveg_shrink[~ np.isnan(wtd_shrink)]
    wtd_1D     = wtd_shrink[~ np.isnan(wtd_shrink)]
    tair_1D    = tair_shrink[~ np.isnan(wtd_shrink)]
    tmax_1D    = tmax_shrink[~ np.isnan(wtd_shrink)]
    tmin_1D    = tmin_shrink[~ np.isnan(wtd_shrink)]
    qle_1D     = qle_shrink[~ np.isnan(wtd_shrink)]
    qh_1D      = qh_shrink[~ np.isnan(wtd_shrink)]
    tveg_1D    = tveg_shrink[~ np.isnan(wtd_shrink)]
    evap_1D    = evap_shrink[~ np.isnan(wtd_shrink)]
    vegt_1D    = vegt_shrink[~ np.isnan(wtd_shrink)]
    fw_1D      = fw_shrink[~ np.isnan(wtd_shrink)]

    df         = pd.DataFrame(iveg_1D, columns=['iveg'])
    df['wtd']  = wtd_1D
    df['tair'] = tair_1D
    df['tmax'] = tmax_1D
    df['tmin'] = tmin_1D
    df['qle']  = qle_1D
    df['qh']   = qh_1D
    df['tveg'] = tveg_1D
    df['evap'] = evap_1D
    df['vegt'] = vegt_1D
    df['fw']   = fw_1D

    # df['rank']     = df[var_name].rank()
    # df_sort        = df.sort_values(by=var_name)

    # ============ Preset of plotting ============
    # cmap     = cm.Set2(np.arange(5))
    pft      = ["BEF","shrub","grass","crop","barren"]
    order    = ["(a)","(b)","(c)","(d)","(e)"]
    iveg_num = [2, 5, 6, 9, 14]
    cmap     = cm.Set2(np.arange(5))
    markers  = ["o","o","o","o","o"]
    mrk_sz   = 4


    # ============ Plotting ============
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[12,8],sharex=True, sharey=True, squeeze=True)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
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

    cnt = 0

    for i in np.arange(2):
        for j in np.arange(2):
            if cnt < 5:


                xxx_value  = fw_val #fw_val #wtd_val
                var_value  = tmax_val
                clr_value  = tair_val #wtd_val #fw_val
                val_min    = 20 #0.
                val_max    = 50  #15.

                ax[i,j].set_xlim(-0.05,1.)
                # ax[i,j].set_xlim(-0.05,0.8)
                ax[i,j].set_ylim(-2.49,0.6)

                if i == 1:
                    ax[i,j].set_xlabel("Δ$β$ (-)", fontsize=12)
                    # ax[i,j].set_xlabel("Δ$β$ (-)", fontsize=12)
                if j == 0:
                    ax[i,j].set_ylabel("ΔTmax ($\mathregular{^o}$C)", fontsize=12)

                # gam fitting
                gam    = LinearGAM(n_splines=4).gridsearch(xxx_value[:,None], var_value) # n_splines=22
                x_pred = np.linspace(min(xxx_value), max(xxx_value), num=100)
                y_pred = gam.predict(x_pred)
                y_int  = gam.confidence_intervals(x_pred, width=.95)

                color_map = get_cmap("gist_heat").reversed()
                sct    = ax[i,j].scatter(xxx_value, var_value, s=mrk_sz, vmin=val_min, vmax=val_max, marker=markers[cnt],c=clr_value, cmap=color_map, label=pft[cnt]) # alpha=0.7,
                ax[i,j].axhline(y=0, color="gray", linestyle='--')
                if cnt == 3:
                    cbar   = plt.colorbar(sct, ax=ax, ticklocation="right", pad=0.01, orientation="vertical", aspect=40, shrink=1.)
                        #  norm=colors.NoNorm(),boundaries=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7], ticks=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7]) # cax=cax,
                    cbar.set_label('T ($\mathregular{^o}$C)',  loc='center') # rotation=270,
                    # cbar.set_label('WTD (m)',  loc='center') # rotation=270,
                ax[i,j].text(0.02, 0.95, order[cnt]+" "+pft[cnt], transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)


                # plt.title(pft[i], size=12)
                # clb.set_clim(0, 1)

                ax[i,j].plot(x_pred, y_pred, color="black", ls='-', lw=2.0, zorder=10)
                ax[i,j].fill_between(x_pred, y_int[:, 0], y_int[:, 1], alpha=0.2,
                                facecolor='grey', zorder=10)

                # ax[i,j].legend()

                cnt = cnt + 1

    fig.savefig("./plots/figures/scatter_delta_Tmax_FWsoil_Tair_PFT_"+message,bbox_inches='tight')#, pad_inches=0.1)

if __name__ == "__main__":

    # =============================== Operation ================================

    # ================================
    # Decks for plot_spatial_land_days
    # ================================
    case_names = ["hw2009_3Nov", "hw2013_3Nov", "hw2019_3Nov"] #
    periods    = ["20090122-20090213","20121229-20130122","20190108-20190130"] #
    time_ss    = [
                   datetime(2009,1,28,0,0,0,0), datetime(2013,1,4,0,0,0,0), datetime(2019,1,14,0,0,0)
                   ] #
    time_es    = [
                   datetime(2009,2,8,23,59,0,0), datetime(2013,1,18,23,59,0,0), datetime(2019,1,26,23,59,0,0)
                   ] #

    message = "GW-FD"

    plot_spatial_land_days(case_names, periods, time_ss, time_es, message=message)


    # # ================================
    # # Decks for plot_spatial_land
    # #           plot_spatial_tair_fwsoil
    # #           plot_spatial_fwsoil_tair
    # # ================================
    # case_names = ["hw2009_3Nov", "hw2013_3Nov", "hw2019_3Nov"]
    # seconds    = None #[18.*60.*60.,6.*60.*60.]
    # for case_name in case_names:

    #     if case_name == "hw2009_3Nov":
    #         period     = "20090122-20090213"
    #         time_s = datetime(2009,1,28,0,0,0,0)
    #         time_e = datetime(2009,2,8,23,59,0,0)
    #     elif  case_name == "hw2013_3Nov":
    #         period     = "20121229-20130122"
    #         time_s = datetime(2013,1,4,0,0,0,0)
    #         time_e = datetime(2013,1,18,23,59,0,0)
    #     elif  case_name == "hw2019_3Nov":
    #         period     = "20190108-20190130"
    #         time_s = datetime(2019,1,14,0,0,0)
    #         time_e = datetime(2019,1,26,23,59,0,0)

    #     cpl_land_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'

    #     cpl_land_file_gw  = cpl_land_file + '/LIS.CABLE.'+period+'_gw.nc'  # land output of wrf-cable run
    #     cpl_land_file_fd  = cpl_land_file + '/LIS.CABLE.'+period+'_fd.nc'  # land output of wrf-cable run

    #     file_paths        = [cpl_land_file_fd,cpl_land_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw

    #     if seconds == None:
    #         if len(file_paths) > 1:
    #             message = "GW-FD_"+str(time_s)+"-"+str(time_e)
    #         else:
    #             message = "GW_"+str(time_s)+"-"+str(time_e)
    #     else:
    #         if seconds[0] < seconds[1]:
    #             day_or_night = "Day"
    #         else:
    #             day_or_night = "Night"
    #         if len(file_paths) > 1:
    #             message = day_or_night+"_GW-FD_"+str(time_s)+"-"+str(time_e)
    #         else:
    #             message = day_or_night+"_Land_GW_"+str(time_s)+"-"+str(time_e)

    #     plot_spatial_land(file_paths, time_s, time_e, seconds=seconds, message=message)
    #     # plot_spatial_fwsoil_tair(file_paths, time_s, time_e, seconds=seconds, message=message)
