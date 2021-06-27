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
sns.set_theme(style="ticks", palette="pastel")

def clean_files(fname, start_yr):

    df = pd.read_csv(fname, header=None)
    df.columns = ["day", "var"]

    # Add correct timestamps
    start = date(start_yr,1,1)
    dates = []
    for i in range(len(df)):
        # need the -1 or the first day is 17 not 16
        delta = timedelta(int(df.day[i]) - 1)
        offset = start + delta
        dates.append(offset)

    df["dates"] = dates
    #df = df.set_index('dates')
    df.index = pd.to_datetime(df.dates, format = '%Y/%m/%d', utc=False)
    df["year"] = df.index.year

    ofname = fname.replace("txt", "csv")
    df.to_csv(ofname, index=True)

def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)

def Fig3_boxplot():
    fname = "EF_FD_rawdata_4_Python.csv"
    df_fd = pd.read_csv(fname)
    df_fd['experiment'] = "FD"
    fname = "EF_GW_rawdata_4_Python.csv"
    df_gw = pd.read_csv(fname)
    df_gw['experiment'] = "GW"

    # you can omit, or do something else with this bit later
    df_fd.temp = np.where(df_fd.temp > 3.8, df_fd.temp == 3.8, df_fd.temp)
    df_fd.temp = np.where(df_fd.temp < 0.0, df_fd.temp == 0.0, df_fd.temp)
    df_gw.temp = np.where(df_gw.temp > 3.8, df_gw.temp == 3.8, df_gw.temp)
    df_gw.temp = np.where(df_gw.temp < 0.0, df_gw.temp == 0.0, df_gw.temp)

    medians_gw = df_gw.groupby(['year'])['temp'].median().values
    medians_fd = df_fd.groupby(['year'])['temp'].median().values


    # make one dataframe
    df = pd.concat([df_fd, df_gw])
    df.drop(['dates.1'], axis=1)

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
    ax2 = ax.twinx()

    ax = sns.boxplot(x="year", y="temp", data=df,
                     showfliers=False, palette=["m", "g"],
                     hue="experiment", whis=0)
    #ax.fill_between(np.arange(2000, 2010), 0, 1, facecolor="lightgrey", alpha=0.5)
    ax2.plot(medians_gw-medians_fd, ls="-", color="black", label="GW-FD")
    #align_yaxis(ax, 0, ax2, 0)

    ax.set_ylim(0, 1.0)
    #ax2.set_ylim(0, 0.2)

    xxlim = ax.get_xlim()
    #print(xxlim)
    ax.fill_between([-0.5, 9], 0, 0.8, facecolor="lightgrey", alpha=0.3)
    ax.fill_between([16.5, 19.5], 0, 0.8, facecolor="lightgrey", alpha=0.3)
    #ax.axhspan(0, 1, facecolor="lightgrey", alpha=0.5)

    #ax.axvspan(2000, 2009, 0, 1, facecolor="lightgrey", alpha=0.5)

    ax.set_ylabel('EF ($\degree$C)')
    ax2.set_ylabel('Temperature ($\degree$C)')
    ax.set_xlabel(" ")
    #ax.legend(numpoints=1, loc="best", frameon=False)
    ax.legend(numpoints=1, loc=(0.7, 0.8), frameon=False)

    for ind, label in enumerate(ax.get_xticklabels()):
        if ind % 2 == 0:  # every 2nd label
            label.set_visible(True)
        else:
            label.set_visible(False)

    of = "/Users/mdekauwe/Desktop/blah.png"
    fig.savefig(of, bbox_inches='tight', dpi=300, pad_inches=0.1)
