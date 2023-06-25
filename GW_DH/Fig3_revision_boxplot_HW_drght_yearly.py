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

    fig, axs = plt.subplots(2, 2, figsize=(12,7), constrained_layout = True)

    # fig = plt.figure(figsize=(12,6))
    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Times New Roman"
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

    orders  = ['(a)','(b)','(c)','(d)']

    for i, var_name in enumerate(var_names):
        print(var_name)

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


        axs2 = axs[row,col].twinx()


        # Adding shadings
        # fill_color = (1., 0.972549, 0.862745) # named color "cornsilk" in ncl
        axs[row,col].fill_between([-0.5, 8.5], ranges[i][0], ranges[i][1], facecolor="gray", alpha=0.1)
        axs[row,col].fill_between([15.5, 18.5], ranges[i][0], ranges[i][1], facecolor="gray", alpha=0.1)

        # Plotting boxplot
        sns.boxplot(x="year", y="var", data=df, ax=axs[row,col],
                        showfliers=False, palette=["b", "r"],
                        hue="experiment", whis=0) # ["m", "g"]

        xtickslocs     = [    0,   1,      2,  3,      4,  5,      6,  7,      8,  9,
                             10,  11,     12,  13,    14, 15,     16, 17,     18     ]
        xticklabels    = ["2001", "", "2003", "", "2005", "", "2007", "", "2009", "",
                          "2011", "", "2013", "", "2015", "", "2017", "", "2019"     ]
        # plt.setp(axs[row,col].get_xticklabels(), visible=False)


        axs[row,col].set_xlim([-1.,19.])
        axs[row,col].set_ylim(ranges[i])
        axs[row,col].set_ylabel(ylabels[i],fontsize=14)
        axs[row,col].yaxis.set_tick_params(labelsize=12)

        axs[row,col].set_xlabel(" ",fontsize=14)
        axs[row,col].set(xticks=xtickslocs, xticklabels=xticklabels)
        axs[row,col].xaxis.set_tick_params(labelsize=12,labelrotation=45)
        # axs[row,col].xtick.label.set_rotation(45)
        # axs[row,col].yaxis.label.set_size(14)
        # axs[row,col].ytick.label.set_size(14)


        if row == 0:
            # axs[row,col].get_xaxis().set_visible(False)
            plt.setp(axs[row,col].get_xticklabels(), visible=False)

        axs[row,col].text(0.01, 0.95, orders[i],transform=axs[row,col].transAxes, fontsize=16, verticalalignment='top', bbox=props) #


        # -------------------- lines ---------------------
        # read line values
        medians_gw = df_gw.groupby(['year'])['var'].median().values
        medians_fd = df_fd.groupby(['year'])['var'].median().values

        # Plotting boxplot
        axs2.plot(medians_gw-medians_fd, ls="-", color="gray", label="GW-FD")
        #align_yaxis(ax, 0, ax2, 0)
        print(medians_gw-medians_fd)
        axs2.set_ylim(ranges_diff[i])
        axs2.set_ylabel(ylabels_R[i],fontsize=14)
        axs2.yaxis.set_tick_params(labelsize=12)
        #
        # axs2.xaxis.label.set_size(14)
        # axs2.xtick.label.set_size(14)
        # set_xticklabels(ax1_x, fontsize=15)
        # axs2.yaxis.label.set_size(14)
        # axs2.ytick.label.set_size(14)

        # for ind, label in enumerate(axs[row,col].get_xticklabels()):
        #     if ind % 2 == 0:  # every 2nd label
        #         label.set_visible(True)
        #         plt.setp(axs[row,col].get_xticklabels(), visible=True)
        #         axs[row,col].set(xticks=xtickslocs, xticklabels=xticklabels)
        #     else:
        #         label.set_visible(False)
        #         plt.setp(axs[row,col].get_xticklabels(), visible=True)
        #         axs[row,col].set(xticks=xtickslocs, xticklabels=xticklabels)

        if row == 0 and col ==0:
            axs[row,col].legend(numpoints=1, loc="best", frameon=False) # loc=(0.7, 0.8)
            # axs2.legend(numpoints=1, loc="best", frameon=False)
        else:
            axs[row,col].get_legend().remove()

    fig.savefig("./plots/Fig3_revision_boxplot_HW_drght_yearly.png", bbox_inches='tight', dpi=300, pad_inches=0.1) #


if __name__ == "__main__":


    var_names   = ["deltaT","EF","TVeg","Fwsoil"]
    ylabels     = ["ΔT ($\mathregular{^o}$C)","EF (-)","Et (mm d$\mathregular{^{-1}}$)","$β$ (-)"]
    ylabels_R   = ["ΔT$\mathregular{_{GW}}$ - ΔT$\mathregular{_{FD}}$ ($\mathregular{^o}$C)","ΔEF (-)","ΔEt (mm d$\mathregular{^{-1}}$)","Δ$β$ (-)"]
    ranges      = [[-0.5, 5], [0., 0.8], [0., 3.8], [0., 1.05]]
    ranges_diff = [[-0.8, 0.], [0., 0.2], [0., 1.1], [0., 0.4]]
    start_yr    = 2000

    Fig3_boxplot(start_yr,var_names,ylabels,ylabels_R,ranges,ranges_diff)
