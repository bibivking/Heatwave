#!/usr/bin/env python

"""
Example plot...

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (24.06.2021)"
__email__ = "mdekauwe@gmail.com"


import pandas as pd
from datetime import date, timedelta
import sys
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def read_summer_heatwave(fname, start_yr):

    df = pd.read_csv(fname, header=None)
    df.columns = ["day", "var"]

    # Add correct timestamps
    start = date(start_yr,1,1)
    dates = []
    for i in range(len(df)):

        delta = timedelta(int(df.day[i]) - 1) # days from start_yr-01-01
        offset = start + delta
        dates.append(offset)

    df["dates"] = dates
    #df = df.set_index('dates')

    df.index = pd.to_datetime(df.dates, format = '%Y/%m/%d', utc=False)
    df["year"] = df.index.year
    df_summer_hw = df[np.any([ np.all( [ df.index.month == 12, df.index.year < 2019 ], axis=0), 
                               np.all( [ df.index.month <= 2,  df.index.year > 2000 ], axis=0) ], axis=0)]
    
    for year in np.arange(2000,2019):
        # df_summer_hw[np.all([df_summer_hw.index.month == 12, df_summer_hw.index.year == year], axis=0)]['year'] = year + 1
        df_summer_hw['year'].loc[np.all([df_summer_hw.index.month == 12, df_summer_hw.index.year == year], axis=0)] = year + 1
       
    print(df_summer_hw)
    return df_summer_hw

def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)

def Fig3_boxplot(start_yr,var_names,ylabels,ylabels_R,ranges,ranges_diff):

    # set plots
    sns.set_style("ticks")
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    sns.set_theme(style="ticks", palette="pastel")

    fig, axs = plt.subplots(2, 2, figsize=(10, 7))

    # fig = plt.figure(figsize=(12,6))
    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 12
    plt.rcParams['font.size']       = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
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

    # ax = fig.add_subplot(111)
    # ax2 = ax.twinx()

    orders  = ['(a)','(b)','(c)','(d)']

    for i, var_name in enumerate(var_names):

        row = i//2 # round
        col = i%2  # mod

        # -------------------- boxplot ---------------------

        # read box values
        filename_GW = "./txt/"+var_name+"_GW_rawdata_4_Python.txt"
        filename_FD = "./txt/"+var_name+"_FD_rawdata_4_Python.txt"

        df_gw = read_summer_heatwave(filename_GW, start_yr)
        df_gw['experiment'] = "GW"
        df_fd = read_summer_heatwave(filename_FD, start_yr)
        df_fd['experiment'] = "FD"

        # make one dataframe
        df = pd.concat([df_gw,df_fd])
        print(df)

        # Plotting boxplot
        axs[row,col] = sns.boxplot(x="year", y="var", data=df,
                        showfliers=False, palette=["m", "g"],
                        hue="experiment", whis=0)
        # xxlim = axs[row,col].get_xlim()
        axs[row,col].set_ylim(ranges[i])
        axs[row,col].set_ylabel(ylabels[i])
        # axs[row,col].set_xlabel(" ")

        # Adding shadings
        fill_color = (1., 0.972549, 0.862745) # named color "cornsilk" in ncl
        axs[row,col].fill_between([ 0.5, 8.5], ranges[i][0], ranges[i][1], facecolor=fill_color, alpha=0.5)
        axs[row,col].fill_between([16.5, 19.5], ranges[i][0], ranges[i][1], facecolor=fill_color, alpha=0.5)

        xtickslocs     = [    0,   1,      2,  3,      4,  5,      6,  7,      8,  9, 
                             10,  11,     12,  13,    14, 15,     16, 17,     18     ]
        xticklabels    = ["2001", "", "2003", "", "2005", "", "2007", "", "2009", "",
                          "2011", "", "2013", "", "2015", "", "2017", "", "2019"     ]
        # plt.setp(axs[row,col].get_xticklabels(), visible=False)
        
        if row == 0 and col ==0:
            axs[row,col].legend(numpoints=1, loc="best", frameon=False) # loc=(0.7, 0.8) 
        else:
            axs[row,col].get_legend().remove()

        if row == 1:    
            # plt.setp(axs[row,col].get_xticklabels(), visible=False)
            # axs[row,col].get_xaxis().set_visible(True)
            axs[row,col].set(xticks=xtickslocs, xticklabels=xticklabels)
        else:
            axs[row,col].get_xaxis().set_visible(False)

        axs[row,col].text(0.05, 0.95, orders[i],transform=axs[row,col].transAxes, fontsize=14, verticalalignment='top', bbox=props) # 


        # # -------------------- lines ---------------------
        # # read line values
        # medians_gw = df_gw.groupby(['year'])['var'].median().values
        # medians_fd = df_fd.groupby(['year'])['var'].median().values

        # # Plotting boxplot
        # axs2 = axs[row,col].twinx()
        # axs2.plot(medians_gw-medians_fd, ls="-", color=almost_black, label="GW-FD")
        # #align_yaxis(ax, 0, ax2, 0)

        # axs2.set_ylim(ranges_diff[i])
        # axs2.set_ylabel(ylabels_R[i])
       
        # for ind, label in enumerate(axs[row,col].get_xticklabels()):
        #     if ind % 2 == 0:  # every 2nd label
        #         label.set_visible(True)
        #         plt.setp(axs[row,col].get_xticklabels(), visible=True)
        #         axs[row,col].set(xticks=xtickslocs, xticklabels=xticklabels)
        #     else:
        #         label.set_visible(False)
        #         plt.setp(axs[row,col].get_xticklabels(), visible=True)
        #         axs[row,col].set(xticks=xtickslocs, xticklabels=xticklabels)

    fig.savefig("./plots/plot_boxplots.png", bbox_inches='tight', dpi=300, pad_inches=2) # 
 

if __name__ == "__main__":


    var_names   = ["deltaT","EF","TVeg","Fwsoil"]
    ylabels     = ["ΔT ($^{o}$C)","EF (-)","Et (mm d-1)","β (-)"]
    ylabels_R   = ["ΔT_GW - ΔT_FD ($^{o}$C)","ΔEF (-)","ΔEt (mm d-1)","Δβ (-)"]
    ranges      = [[-0.5, 5], [0., 0.8], [0., 3.8], [0., 1.05]]
    ranges_diff = [[-0.8, 0.], [0., 0.2], [0., 1.1], [0., 0.4]]
    start_yr    = 2000

    Fig3_boxplot(start_yr,var_names,ylabels,ylabels_R,ranges,ranges_diff)
