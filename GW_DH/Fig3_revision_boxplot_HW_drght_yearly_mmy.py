#!/usr/bin/env python

__author__  = "MU Mengyuan"
__version__ = "1.0 (22.03.2021)"
__email__   = "mu.mengyuan815@gmail.com"

import os
import sys
import glob
import numpy as np
from numpy.core.numeric import NaN
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import datetime as dt
import netCDF4 as nc
import scipy.stats as stats
import seaborn as sns
from matplotlib import cm
from matplotlib import ticker
from scipy.interpolate import griddata
from sklearn.metrics import mean_squared_error

def Fig3_boxplot(var_names,ylabels,ylabels_R,ranges,ranges_diff):

    """

    """

    # ======================= Plot setting ============================
    fig, axs = plt.subplots(2, 2, figsize=(7, 5))
    fig.subplots_adjust(hspace=0.20)
    fig.subplots_adjust(wspace=0.12)
    print(axs)
    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams["legend.markerscale"] = 3.0

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']  = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor']  = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    # set the box type of sequence number
    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
    #colors = cm.Set2(np.arange(0,len(case_labels)))

    # ==================== Summers ===================
            # 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,\
            # 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019
    ts_s    = [ 335, 700, 1065, 1430, 1796, 2161, 2526, 2891, 3257, 3622,
                3987,4352, 4718, 5083, 5448, 5813, 6179, 6544, 6909]

    ts_e    = [ 425, 790, 1155, 1521, 1886, 2251, 2616, 2982, 3347, 3712,
                4077, 4443, 4808,5173, 5538, 5904, 6269, 6634, 6999]

    orders  = ['(a)','(b)','(c)','(d)']

    # ==================== Start ====================
    for i, var_name in enumerate(var_names):

        row = i//2 # round
        col = i%2  # mod

        filename_GW = "./txt/"+var_name+"_GW_rawdata_4_Python.txt"
        filename_FD = "./txt/"+var_name+"_FD_rawdata_4_Python.txt"

        var_GW      = pd.read_csv(filename_GW, header = None, names= ["date","GW"])
        var_FD      = pd.read_csv(filename_FD, header = None, names= ["date","FD"])
        print(var_GW)
        data_len    = len(var_GW)
        var_all     = pd.DataFrame(np.zeros(data_len*2),columns=['values'])

        var_all['values'][0:data_len]          = var_GW["GW"]
        var_all['values'][data_len:data_len*2] = var_FD["FD"]

        var_all['case']                        = [''] * data_len*2
        var_all['case'][0:data_len]            = ['GW'] * data_len
        var_all['case'][data_len:data_len*2]   = ['FD'] * data_len

        var_all['date']                        = [0] * data_len*2
        var_all['date'][0:data_len]            = var_GW["date"]
        var_all['date'][data_len:data_len*2]   = var_FD["date"]

        var_all['year']                        = NaN * data_len*2

        for k in np.arange(len(var_all)):
            for j in np.arange(len(ts_s)):
                if var_all['date'][k] >= ts_s[j] and var_all['date'][k] <= ts_e[j]:
                    var_all['year'][k] = 2001+k

        print(var_all['year'])


        # ========================= box-whisker ============================
        # seaborn
        #sns.color_palette("Set2", 8)
        #flatui = ["orange", "dodgerblue"]
        #aaa = sns.set_palette(flatui)
        sns.boxplot(x="year", y="values", hue="case", data=var_all, palette="BrBG",
                    order=np.arange(2001,2019),  width=0.7, hue_order=['GW','FD'],
                    ax=axs[i], showfliers=False, whis=0, color=almost_black) # palette="Set2",

        axs[row,col].set_ylabel(ylabels[i])
        axs[row,col].set_xlabel("")
        axs[row,col].axis('tight')
        axs[row,col].set_ylim(ranges[i])
        axs[row,col].legend(loc='best', frameon=False)

        #colors = cm.Set3(np.arange(0,len(case_labels)))


        # ========================= difference lines ============================
        file_metrics = "./txt/"+var_name+"_CTL_FD_yearly_box_stats.txt"

        var_median   = pd.read_fwf(file_metrics, header = None,
                                   names= ['min','x25','median','x75','max'])
        median_diff  = np.zeros(19)
        print(var_median)

        for l in np.arange(0,19):
            loc1 = l*5+2 # location for GW tree
            loc2 = l*5+3 # location for FD tree

            median_diff[l] = var_median['median'][loc1] - var_median['median'][loc2]

        print(median_diff)

        y = np.arange(1,20)
        print(axs[row,col])
        axs[row,col].plot(y, median_diff,alpha=0.45, c=almost_black)

        axs[row,col].set_xlim()
        axs[row,col].set_ylim(ranges_diff[i])
        axs[row,col].set_ylabel(ylabels_R[i])

        # ========================= shadings ============================
        color = (1., 0.972549, 0.862745) # named color "cornsilk" in ncl
        axs[row,col].axvspan(1,  9, facecolor=color, alpha=0.5)
        axs[row,col].axvspan(17, 19, facecolor=color, alpha=0.5)

        # ========================== order ========================
        axs[row,col].text(0.02, 0.95, orders[i], fontsize=14, verticalalignment='top', bbox=props) # , transform=axs[i].transAxes
        #plt.setp(ax2.get_yticklabels(), visible=False)

    fig.savefig("./plots/Fig3_revision_boxplot" , bbox_inches='tight', pad_inches=0.1)


if __name__ == "__main__":


    var_names   = ["EF"]      #["deltaT","EF","TVeg","Fwsoil"]
    ylabels     = ["EF (-)"]  #["ΔT (oC)","EF (-)","Et (mm d-1)","β (-)"]
    ylabels_R   = ["ΔEF (-)"] #["ΔT_GW - ΔT_FD (oC)","ΔEF (-)","ΔEt (mm d-1)","Δβ (-)"]
    ranges      = [0., 3.8] #[[-0.5, 5], [0., 0.8], [0., 3.8], [0., 1.05]]
    ranges_diff = [0., 1.1] #[[-0.8, 0.], [0., 0.2], [0., 1.1], [0., 0.4]]

    Fig3_boxplot(var_names,ylabels,ylabels_R,ranges,ranges_diff)
