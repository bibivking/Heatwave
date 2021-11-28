#!/usr/bin/python

'''
Plot spitial map of land diagnosis and parameters from LIS-CABLE
1. per time step
2. time period average
'''

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from pygam import LinearGAM
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_range_diff
from common_utils import *

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
                    ax[i,j].set_ylabel("ΔTmax ($\mathregular{^o}$C)", fontsize=12)



                # gam fitting
                gam    = LinearGAM(n_splines=4).gridsearch(xxx_value[:,None], var_value) # n_splines=22
                x_pred = np.linspace(min(xxx_value), max(xxx_value), num=100)
                y_pred = gam.predict(x_pred)
                y_int  = gam.confidence_intervals(x_pred, width=.95)

                sct    = ax[i,j].scatter(xxx_value, var_value, s=mrk_sz, vmin=val_min, vmax=val_max, marker=markers[cnt],c=clr_value, cmap="BrBG", label=pft[cnt]) # alpha=0.7,
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

    fig.savefig("./plots/figures/scatter_delta_Tmax_WTD_PFT_"+message,bbox_inches='tight')#, pad_inches=0.1)

# def plot_spatial_land(file_paths, time_s, time_e, seconds=None,loc_lat=None, loc_lon=None, message=None):
#
#     # Open the NetCDF4 file (add a directory path if necessary) for reading:
#     file1  = Dataset(file_paths[0], mode='r')
#
#     Time   = nc.num2date(file1.variables['time'][:],file1.variables['time'].units,
#             only_use_cftime_datetimes=False, only_use_python_datetimes=True)
#     time   = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)
#
#     Tair1  = file1.variables['Tair_f_inst']
#     tair1  = spital_var(time,Tair1,time_s,time_e,seconds)
#     tmax1  = spital_var_max(time,Tair1,time_s,time_e,seconds)
#     tmin1  = spital_var_min(time,Tair1,time_s,time_e,seconds)
#
#     Qle1   = file1.variables['Qle_tavg']
#     qle1   = spital_var(time,Qle1,time_s,time_e,seconds)
#     Qh1    = file1.variables['Qh_tavg']
#     qh1    = spital_var(time,Qh1,time_s,time_e,seconds)
#     TVeg1  = file1.variables['TVeg_tavg']
#     tveg1  = spital_var(time,TVeg1,time_s,time_e,seconds)
#     Evap1  = file1.variables['Evap_tavg']
#     evap1  = spital_var(time,Evap1,time_s,time_e,seconds)
#     VegT1  = file1.variables['VegT_tavg']
#     vegt1  = spital_var(time,VegT1,time_s,time_e,seconds)
#     FW1    = file1.variables['FWsoil_tavg']
#     fw1    = spital_var(time,FW1,time_s,time_e,seconds)
#
#     if len(file_paths) > 1:
#         file2  = Dataset(file_paths[1], mode='r')
#         Tair2  = file2.variables['Tair_f_inst']
#
#         tair2  = spital_var(time,Tair2,time_s,time_e,seconds)
#         tmax2  = spital_var_max(time,Tair2,time_s,time_e,seconds)
#         tmin2  = spital_var_min(time,Tair2,time_s,time_e,seconds)
#
#         Qle2   = file2.variables['Qle_tavg']
#         qle2   = spital_var(time,Qle2,time_s,time_e,seconds)
#         Qh2    = file2.variables['Qh_tavg']
#         qh2    = spital_var(time,Qh2,time_s,time_e,seconds)
#         TVeg2  = file2.variables['TVeg_tavg']
#         tveg2  = spital_var(time,TVeg2,time_s,time_e,seconds)
#         Evap2  = file2.variables['Evap_tavg']
#         evap2  = spital_var(time,Evap2,time_s,time_e,seconds)
#         VegT2  = file2.variables['VegT_tavg']
#         vegt2  = spital_var(time,VegT2,time_s,time_e,seconds)
#         FW2    = file2.variables['FWsoil_tavg']
#         fw2    = spital_var(time,FW2,time_s,time_e,seconds)
#
#         iveg    = file2.variables["Landcover_inst"][0,:,:]
#         WTD_tmp = file2.variables["WaterTableD_tavg"]
#         wtd     = spital_var(time,WTD_tmp,time_s,time_e,seconds)
#
#         tair     = tair2 - tair1
#         tmax     = tmax2 - tmax1
#         tmin     = tmin2 - tmin1
#         qle      = qle2 - qle1
#         qh       = qh2 - qh1
#         tveg     = tveg2 - tveg1
#         evap     = evap2 - evap1
#         vegt     = vegt2 - vegt1
#         fw       = fw2 - fw1
#     else:
#         tair     = tair1
#         tmax     = tmax1
#         tmin     = tmin1
#         qle      = qle1
#         qh       = qh1
#         tveg     = tveg1
#         evap     = evap1
#         vegt     = vegt1
#         fw       = fw1
#
#
#     # # remove ocean and boundary pixels
#     # nlon = len(tair[0,6:-5])
#     # nlat = len(tair[6:-5,0])
#
#     # the new array doesn't set a missing value, and the missing value from old array will become nan
#     # iveg_shrink    = np.zeros(nlat * nlon)
#
#     iveg_shrink    = np.reshape(iveg[6:-5,6:-5],-1)
#     wtd_shrink     = np.reshape(wtd[6:-5,6:-5],-1)
#     tair_shrink    = np.reshape(tair[6:-5,6:-5],-1)
#     tmax_shrink    = np.reshape(tmax[6:-5,6:-5],-1)
#     tmin_shrink    = np.reshape(tmin[6:-5,6:-5],-1)
#     qle_shrink     = np.reshape(qle[6:-5,6:-5],-1)
#     qh_shrink      = np.reshape(qh[6:-5,6:-5],-1)
#     tveg_shrink    = np.reshape(tveg[6:-5,6:-5],-1)
#     evap_shrink    = np.reshape(evap[6:-5,6:-5],-1)
#     vegt_shrink    = np.reshape(vegt[6:-5,6:-5],-1)
#     fw_shrink      = np.reshape(fw[6:-5,6:-5],-1)
#
#     # print("wtd_shrink",~ np.isnan(wtd_shrink))
#     iveg_1D    = iveg_shrink[~ np.isnan(wtd_shrink)]
#     wtd_1D     = wtd_shrink[~ np.isnan(wtd_shrink)]
#     tair_1D    = tair_shrink[~ np.isnan(wtd_shrink)]
#     tmax_1D    = tmax_shrink[~ np.isnan(wtd_shrink)]
#     tmin_1D    = tmin_shrink[~ np.isnan(wtd_shrink)]
#     qle_1D     = qle_shrink[~ np.isnan(wtd_shrink)]
#     qh_1D      = qh_shrink[~ np.isnan(wtd_shrink)]
#     tveg_1D    = tveg_shrink[~ np.isnan(wtd_shrink)]
#     evap_1D    = evap_shrink[~ np.isnan(wtd_shrink)]
#     vegt_1D    = vegt_shrink[~ np.isnan(wtd_shrink)]
#     fw_1D      = fw_shrink[~ np.isnan(wtd_shrink)]
#
#     df         = pd.DataFrame(iveg_1D, columns=['iveg'])
#     df['wtd']  = wtd_1D
#     df['tair'] = tair_1D
#     df['tmax'] = tmax_1D
#     df['tmin'] = tmin_1D
#     df['qle']  = qle_1D
#     df['qh']   = qh_1D
#     df['tveg'] = tveg_1D
#     df['evap'] = evap_1D
#     df['vegt'] = vegt_1D
#     df['fw']   = fw_1D
#
#     # df['rank']     = df[var_name].rank()
#     # df_sort        = df.sort_values(by=var_name)
#
#     # ============ Preset of plotting ============
#     # cmap     = cm.Set2(np.arange(5))
#     pft      = ["BEF","shrub","grass","crop","barren"]
#     order    = ["(a)","(b)","(c)","(d)","(e)"]
#     iveg_num = [2, 5, 6, 9, 14]
#     cmap     = cm.Set2(np.arange(5))
#     markers  = ["o","o","o","o","o"]
#     mrk_sz   = 4
#
#
#     # ============ Plotting ============
#     fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[15,10])
#     plt.subplots_adjust(wspace=0.0, hspace=0.0)
#
#     plt.rcParams['text.usetex']     = False
#     plt.rcParams['font.family']     = "sans-serif"
#     plt.rcParams['font.serif']      = "Helvetica"
#     plt.rcParams['axes.linewidth']  = 1.5
#     plt.rcParams['axes.labelsize']  = 14
#     plt.rcParams['font.size']       = 14
#     plt.rcParams['legend.fontsize'] = 12
#     plt.rcParams['xtick.labelsize'] = 12
#     plt.rcParams['ytick.labelsize'] = 14
#
#     almost_black = '#262626'
#     # change the tick colors also to the almost black
#     plt.rcParams['ytick.color']     = almost_black
#     plt.rcParams['xtick.color']     = almost_black
#
#     # change the text colors also to the almost black
#     plt.rcParams['text.color']      = almost_black
#
#     # Change the default axis colors from black to a slightly lighter black,
#     # and a little thinner (0.5 instead of 1)
#     plt.rcParams['axes.edgecolor']  = almost_black
#     plt.rcParams['axes.labelcolor'] = almost_black
#
#     props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
#
#     cnt = 0
#
#     # plot 1
#     i, j = 0, 0
#     for cnt in np.arange(4):
#
#         wtd_val    = df['wtd'][df.iveg == iveg_num[cnt]].values/1000.
#         tair_val   = df['tair'][df.iveg == iveg_num[cnt]].values
#         tmax_val   = df['tmax'][df.iveg == iveg_num[cnt]].values
#         tmin_val   = df['tmin'][df.iveg == iveg_num[cnt]].values
#         qle_val    = df['qle'][df.iveg == iveg_num[cnt]].values
#         qh_val     = df['qh'][df.iveg == iveg_num[cnt]].values
#         tveg_val   = df['tveg'][df.iveg == iveg_num[cnt]].values*3600*24
#         evap_val   = df['evap'][df.iveg == iveg_num[cnt]].values*3600*24
#         vegt_val   = df['vegt'][df.iveg == iveg_num[cnt]].values
#         fw_val     = df['fw'][df.iveg == iveg_num[cnt]].values
#
#         xxx_value  = wtd_val #fw_val #wtd_val
#         var_value  = tmax_val
#         clr_value  = fw_val #wtd_val #fw_val
#         val_min    = -0.1 #0.
#         val_max    = 0.8  #15.
#
#         # gam fitting
#         gam    = LinearGAM(n_splines=4).gridsearch(xxx_value[:,None], var_value) # n_splines=22
#         x_pred = np.linspace(min(xxx_value), max(xxx_value), num=100)
#         y_pred = gam.predict(x_pred)
#         y_int  = gam.confidence_intervals(x_pred, width=.95)
#
#         sct    = ax[i,j].scatter(xxx_value, var_value, s=mrk_sz, vmin=val_min, vmax=val_max, marker=markers[cnt],c=clr_value, cmap="BrBG", label=pft[cnt]) # alpha=0.7,
#         cbar   = plt.colorbar(sct, ax=ax[i,j], ticklocation="right", pad=0.01, orientation="vertical", aspect=20, shrink=0.6)
#                 #  norm=colors.NoNorm(),boundaries=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7], ticks=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7]) # cax=cax,
#
#         cbar.set_label('Δ$β$ (-)',  loc='center') # rotation=270,
#         # cbar.set_label('WTD (m)',  loc='center') # rotation=270,
#         ax[i,j].text(0.02, 0.95, order[cnt]+" "+pft[cnt], transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)
#
#         # plt.title(pft[i], size=12)
#         # clb.set_clim(0, 1)
#
#         ax[i,j].plot(x_pred, y_pred, color="black", ls='-', lw=2.0, zorder=10)
#         ax[i,j].fill_between(x_pred, y_int[:, 0], y_int[:, 1], alpha=0.2,
#                         facecolor='grey', zorder=10)
#
#         # ax[i,j].legend()
#         ax[i,j].set_xlabel("WTD (m)", fontsize=12)
#         # ax[i,j].set_xlabel("Δ$β$ (-)", fontsize=12)
#         ax[i,j].set_ylabel("ΔTmax ($\mathregular{^o}$C)", fontsize=12)
#
#         ax[i,j].set_xlim(-0.,15)
#         # ax[i,j].set_xlim(-0.05,0.8)
#         ax[i,j].set_ylim(-2.8,0.5)
#
#
#     fig.savefig("./plots/figures/scatter_delta_Tmax_WTD_PFT_"+message,bbox_inches='tight')#, pad_inches=0.1)

if __name__ == "__main__":


    # =============================== Operation ================================
    case_names = ["hw2009_3Nov", "hw2013_3Nov", "hw2019_3Nov"]
    seconds    = None #[18.*60.*60.,6.*60.*60.]
    for case_name in case_names:

        if case_name == "hw2009_3Nov":
            period     = "20090122-20090213"
            time_s = datetime(2009,1,28,0,0,0,0)
            # time_e = datetime(2009,1,28,23,59,0,0)
            time_e = datetime(2009,2,8,23,59,0,0)
        elif  case_name == "hw2013_3Nov":
            period     = "20121229-20130122"
            time_s = datetime(2013,1,4,0,0,0,0)
            time_e = datetime(2013,1,18,23,59,0,0)
        elif  case_name == "hw2019_3Nov":
            period     = "20190108-20190130"
            time_s = datetime(2019,1,14,0,0,0)
            time_e = datetime(2019,1,26,23,59,0,0)

        cpl_land_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'

        cpl_land_file_gw  = cpl_land_file + '/LIS.CABLE.'+period+'_gw.nc'  # land output of wrf-cable run
        cpl_land_file_fd  = cpl_land_file + '/LIS.CABLE.'+period+'_fd.nc'  # land output of wrf-cable run

        file_paths        = [cpl_land_file_fd,cpl_land_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw

        if seconds == None:
            if len(file_paths) > 1:
                message = "GW-FD_"+str(time_s)+"-"+str(time_e)
            else:
                message = "GW_"+str(time_s)+"-"+str(time_e)
        else:
            if seconds[0] < seconds[1]:
                day_or_night = "Day"
            else:
                day_or_night = "Night"
            if len(file_paths) > 1:
                message = day_or_night+"_GW-FD_"+str(time_s)+"-"+str(time_e)
            else:
                message = day_or_night+"_Land_GW_"+str(time_s)+"-"+str(time_e)

        plot_spatial_land(file_paths, time_s, time_e, seconds=seconds, message=message)
