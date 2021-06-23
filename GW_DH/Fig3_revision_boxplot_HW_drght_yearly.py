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
    (a) box-whisker of fwsoil
    (b) fwsoil vs SM
    """

    # ======================= Plot setting ============================
    fig, axs = plt.subplots(2, 2, figsize=(10, 11))
    fig.subplots_adjust(hspace=0.20)
    fig.subplots_adjust(wspace=0.12)

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

        filename_GW = "./txt/"+var_name+"_GW_rawdata_4_Python.txt"
        filename_FD = "./txt/"+var_name+"_FD_rawdata_4_Python.txt"

        var_GW      = pd.read_csv(filename_GW, header = None, names= ["date","GW"])
        var_FD      = pd.read_csv(filename_FD, header = None, names= ["date","FD"])

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

        year = 2001
        for i in np.arange(len(ts_s)):
            var_all['year'][var_all['date'] >= ts_s[i] and var_all['date'] <= ts_e[i]] = year+i


        # ========================= box-whisker ============================
       
        # seaborn
        #sns.color_palette("Set2", 8)
        #flatui = ["orange", "dodgerblue"]
        #aaa = sns.set_palette(flatui)
        sns.boxplot(x="year", y="values", hue="case", data=var_all, palette="BrBG",
                    order=np.arange(2001,2019),  width=0.7, hue_order=['GW','FD'],
                    ax=axs[i], showfliers=False, color=almost_black) # palette="Set2",

        axs[i].set_ylabel(ylabels[i])
        # ax[i].set_ylabel(ylabels[i])

        axs[i].set_xlabel("")
        axs[i].axis('tight')
        axs[i].set_ylim(ranges[i])
        axs[i].legend(loc='best', frameon=False)
        axs[i].text(0.02, 0.95, , transform=axs[i].transAxes, fontsize=14, verticalalignment='top', bbox=props)

        #colors = cm.Set3(np.arange(0,len(case_labels)))



        # ========================= difference lines ============================
        color = "black"
        y = np.arange(1,19)
        diff = ???
        ax[i].scatter(y, diff,  s=1., marker='-', alpha=0.45, c=color)

        ax[i].set_xlim()
        ax[i].set_ylim(ranges_diff[i])
        ax[i].set_ylabel(ylabels_R[i])

        ax[i].legend(numpoints=1, loc='best', frameon=False)
        ax[i].text(0.02, 0.95, orders[i], transform=ax[i].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        #plt.setp(ax2.get_yticklabels(), visible=False)

    fig.savefig("./plots/Fig3_revision_boxplot" , bbox_inches='tight', pad_inches=0.1)



'''
Martin's script 

#!/usr/bin/env python

"""
Plot SWP

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (18.10.2017)"
__email__ = "mdekauwe@gmail.com"

import xarray as xr
import matplotlib.pyplot as plt
import sys
import datetime as dt
import pandas as pd
import numpy as np
from matplotlib.ticker import FixedLocator
import datetime
import os
import glob
import seaborn as sns
sns.set_theme(style="ticks", palette="pastel")


def main(fname, fname2, fname_fixed):

    df_fix = pd.read_csv(fname_fixed)
    df_fix = df_fix[df_fix.plc_max > -500.].reset_index()
    df_fix['map'] = np.ones(len(df_fix)) * np.nan

    df = pd.read_csv(fname)
    df = df[df.plc_max > -500.].reset_index()
    df['map'] = np.ones(len(df)) * np.nan

    df_map = pd.read_csv("euc_map.csv")
    df_map = df_map.sort_values(by=['map'])
    print(df_map)

    species = df_map.species
    species = species.str.replace("_", " ")
    species = species.str.replace("Eucalyptus", "E.")


    df2 = pd.read_csv(fname2)
    df2['map'] = np.ones(len(df2)) * np.nan

    for i in range(len(df_map)):
        spp = df_map.species[i]
        map = df_map.map[i]

        for j in range(len(df)):

            if df.species[j] == spp.replace(" ", "_"):
                df['map'][j] = map
                df_fix['map'][j] = map

    for i in range(len(df_map)):
        spp = df_map.species[i]
        map = df_map.map[i]

        for k in range(len(df2)):
            if df2.species[k] == spp:
                df2['map'][k] = map
    df2 = df2.sort_values(by=['map'])
    sorted_map = df_map['map'].values

    #df = df.sort_values(by=['map'])

    #for i in range(len(df)):
    #    print(i, df.species[i], df.psi_leaf_mean[i])

    df['species'] = df['species'].str.replace("_", " ")
    df['species'] = df['species'].str.replace("Eucalyptus", "E.")

    df_fix['species'] = df_fix['species'].str.replace("_", " ")
    df_fix['species'] = df_fix['species'].str.replace("Eucalyptus", "E.")

    df2['species'] = df2['species'].str.replace("_", " ")
    df2['species'] = df2['species'].str.replace("Eucalyptus", "E.")
    #species = np.unique(df.species)


    sns.set_style("ticks")
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    fig = plt.figure(figsize=(12,6))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax = fig.add_subplot(111)

    flierprops = dict(marker='o', markersize=3, markerfacecolor="black")


    #ax = sns.boxplot(x="species", y="psi_leaf_min", data=df,
    #                 flierprops=flierprops, palette=["m", "g", "b"],
    #                 hue="time", order=species)
    lineprops = {'color': 'black', 'linewidth': 2.0}
    boxplot_kwargs = dict({'medianprops': lineprops})

    ax = sns.boxplot(x="species", y="psi_leaf_min", data=df,
                     showfliers=False, palette=["m", "g", "b"],
                     hue="time", order=species, **boxplot_kwargs)
    bx = sns.stripplot(x="species", y="psi_leaf_min", data=df,
                       color=".3", size=1.3, dodge=True,
                       hue="time", order=species)

    # fixed co2
    boxprops = {'edgecolor': 'grey', 'linewidth': 1.5, 'facecolor': 'w', 'alpha':.5}
    lineprops = {'color': 'grey', 'linewidth': 1.5}
    whiskerprops = {'ls': ' ', 'lw': 0.0, 'color': 'white'}
    cap_col = dict(color="white")
    medianprops = {'linestyle':'-', 'linewidth':1.5, 'color':'salmon'}
    boxplot_kwargs = dict({'boxprops': boxprops, 'medianprops': medianprops,
                       'whiskerprops': whiskerprops, 'capprops':cap_col})


    cx = sns.boxplot(x="species", y="psi_leaf_min", data=df_fix,
                     showfliers=False, palette=["m", "g", "b"],
                     hue="time", order=species, **boxplot_kwargs)

    ax = sns.scatterplot(data=df2, x="species", y="p50", color="red",
                         marker="*", s=100, label="p50")
    ax.set_ylabel("$\Psi$$_{l}$ (MPa)")
    ax.set_xlabel(" ")

    ##8da0cb
    plt.text(-1.4, -0.2,
             "MAP\n(mm yr$^{-1}$)", horizontalalignment='center', size=10,
             color="royalblue", weight="bold")

    offset = 0
    for i,val in enumerate(sorted_map):

        plt.text(-0.05+offset, -0.2,
                 "%d" % (val), horizontalalignment='center', size=10,
                 color="royalblue")

        offset += 1.0

    #ax.legend(numpoints=1, loc="best", frameon=False)





    handles, labels = ax.get_legend_handles_labels()
    print(handles)
    print(labels)
    hh = handles[0:3]
    ll = labels[0:3]
    hh2 = handles[9:]
    ll2 = labels[9:]

    l = plt.legend(hh+hh2, ll+ll2, loc="best", frameon=False)


    ax.set_xticklabels(species, rotation=90)
    of = "/Users/mdekauwe/Desktop/psi_leaf_all_years.png"
    plt.savefig(of, bbox_inches='tight', dpi=300, pad_inches=0.1)
'''

if __name__ == "__main__":


    var_names   = ["deltaT","EF","TVeg","Fwsoil"]
    ylabels     = ["ΔT (oC)","EF (-)","Et (mm d-1)","β (-)"]
    ylabels_R   = ["ΔT_GW - ΔT_FD (oC)","ΔEF (-)","ΔEt (mm d-1)","Δβ (-)"]
    ranges      = [[-0.5, 5], [0., 0.8], [0., 3.8], [0., 1.05]]
    ranges_diff = [[-0.8, 0.], [0., 0.2], [0., 1.1], [0., 0.4]]
  
    Fig3_boxplot(var_names,ylabels,ylabels_R,ranges,ranges_diff)
