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
from pygam import LinearGAM
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_range_diff
from common_utils import *

def plot_spatial_land(file_paths, var_name, clr_var_name, time_s, time_e, seconds=None,loc_lat=None, loc_lon=None, message=None):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    file1 = Dataset(file_paths[0], mode='r')

    Time  = nc.num2date(file1.variables['time'][:],file1.variables['time'].units,
            only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time  = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)

    Var1  = file1.variables[var_name]
    var1  = spital_var_max(time,Var1,time_s,time_e,seconds)
    
    CLR1  = file1.variables[clr_var_name]
    clr1  = spital_var_max(time,CLR1,time_s,time_e,seconds)

    if len(file_paths) > 1:
        file2   = Dataset(file_paths[1], mode='r')
        Var2    = file2.variables[var_name]
        var2    = spital_var_max(time,Var2,time_s,time_e,seconds)
        
        iveg    = file2.variables["Landcover_inst"][0,:,:]
        WTD_tmp = file2.variables["WaterTableD_tavg"]
        wtd     = spital_var(time,WTD_tmp,time_s,time_e,seconds)

        CLR2    = file2.variables[clr_var_name]
        clr2    = spital_var_max(time,CLR2,time_s,time_e,seconds)

        var     = var2 - var1
        clr     = clr2 - clr1
    else:
        var     = var1

    # remove ocean and boundary pixels
    nlon = len(var[0,6:-5])
    nlat = len(var[6:-5,0])

    # the new array doesn't set a missing value, and the missing value from old array will become nan
    iveg_shrink    = np.zeros(nlat * nlon)
    wtd_shrink     = np.zeros(nlat * nlon)
    var_shrink     = np.zeros(nlat * nlon)
    tmax_shrink    = np.zeros(nlat * nlon)
    clr_shrink      = np.zeros(nlat * nlon)

    iveg_shrink[:] = np.reshape(iveg[6:-5,6:-5],-1)
    wtd_shrink[:]  = np.reshape(wtd[6:-5,6:-5],-1)
    var_shrink[:]  = np.reshape(var[6:-5,6:-5],-1)
    tmax_shrink[:] = np.reshape(var1[6:-5,6:-5],-1)
    clr_shrink[:]  = np.reshape(clr[6:-5,6:-5],-1)

    # print("wtd_shrink",~ np.isnan(wtd_shrink))
    iveg_1D        = iveg_shrink[~ np.isnan(wtd_shrink)]
    wtd_1D         = wtd_shrink[~ np.isnan(wtd_shrink)]
    var_1D         = var_shrink[~ np.isnan(wtd_shrink)]
    tmax_1D        = tmax_shrink[~ np.isnan(wtd_shrink)]
    clr_1D          = clr_shrink[~ np.isnan(wtd_shrink)]

    df             = pd.DataFrame(var_1D, columns=[var_name])
    df['iveg']     = iveg_1D
    df['wtd']      = wtd_1D
    df['tmax']     = tmax_1D
    df['clr']      = clr_1D
    df['rank']     = df[var_name].rank()
    df_sort        = df.sort_values(by=var_name)

    # ==== Plotting ====
    fig     = plt.figure(figsize=(12,9))
    ax      = fig.add_subplot(1,1,1)

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

    # ax.set_xlim(-0.,15)
    ax.set_ylim(-9.,1.)

    iveg_num = [2, 5, 6, 9, 14]
    # cmap     = cm.Set2(np.arange(5))
    cmap     = cm.Set2(np.arange(5))
    pft      = ["BEF","shrub","grass","crop","barren"]
    markers  = ["o","o","o","o","o"]
    mrk_sz   = 4

    i = 4
    # data for plotting
    clr_val     = df_sort['wtd'][df_sort.iveg == iveg_num[i]].values/1000.
                # df_sort['clr'][df_sort.iveg == iveg_num[i]]
                # df_sort['tmax'][df_sort.iveg == iveg_num[i]]-273.15
    clr_val_max = np.nanmax(clr_val)
    clr_val_pct = (clr_val)#/clr_val_max

    xxx_value   = df_sort['clr'][df_sort.iveg == iveg_num[i]]
    # df_sort['wtd'][df_sort.iveg == iveg_num[i]].values/1000.
    var_value   = df_sort[var_name][df_sort.iveg == iveg_num[i]].values

    # gam fitting
    gam    = LinearGAM(n_splines=4).gridsearch(xxx_value[:,None], var_value) # n_splines=22
    x_pred = np.linspace(min(xxx_value), max(xxx_value), num=100)
    y_pred = gam.predict(x_pred)
    y_int  = gam.confidence_intervals(x_pred, width=.95)

    sct = ax.scatter( xxx_value, var_value, s=mrk_sz, marker=markers[i], 
                c=clr_val_pct, cmap="Blues", label=pft[i]) # alpha=0.7, 
    plt.colorbar(sct)

    ax.plot(x_pred, y_pred, color="black", ls='-', lw=2.0, zorder=10)
    ax.fill_between(x_pred, y_int[:, 0], y_int[:, 1], alpha=0.2,
                    facecolor='grey', zorder=10)
    ax.legend()

    # cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    # ax.text(0.03, 0.95, '(a)', transform=ax.transAxes, fontsize=18, verticalalignment='top', bbox=props)                
    

    fig.savefig("./plots/5Nov/scatter_deltaT_WTD_PFT-"+pft[i],bbox_inches='tight')#, pad_inches=0.1)

if __name__ == "__main__":


    # =============================== Operation ================================
    case_name  = "hw2009_3Nov" # "hw2013_3Nov"# "hw2019_3Nov"#

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
    seconds           = None #[18.*60.*60.,6.*60.*60.]

    var_name = 'VegT_tavg' #'Tair_f_inst' #"Qle_tavg" #,"Qle_tavg"#'Tair_f_inst'
    clr_var_name = 'FWsoil_tavg'

    if len(file_paths) > 1:
        message = "Land_GW-FD_"+str(time_s)+"-"+str(time_e)
    else:
        message = "Land_GW_"+str(time_s)+"-"+str(time_e)

    plot_spatial_land(file_paths, var_name, clr_var_name, time_s, time_e, seconds=seconds, message=message)
