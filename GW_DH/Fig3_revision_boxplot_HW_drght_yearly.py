#!/usr/bin/env python

__author__  = "MU Mengyuan"
__version__ = "1.0 (22.03.2021)"
__email__   = "mu.mengyuan815@gmail.com"

import os
import sys
import glob
import numpy as np
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
from plot_eucface_get_var import *

def plot_fwsoil_boxplot_SM( fcables, case_labels, layers, ring):

    """
    (a) box-whisker of fwsoil
    (b) fwsoil vs SM
    """

    # ======================= Plot setting ============================
    fig = plt.figure(figsize=[10,11])
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

    # choose colormap
    #colors = cm.Set2(np.arange(0,len(case_labels)))

    #print(colors)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    # ========================= box-whisker of fwsoil============================
    day_start_drought = 2101 # 2017-10-1
    day_end_drought   = 2466 # 2018-10-1

    day_start_all     = 367  # first day of 2013
    day_end           = 2923 # first day of 2020

    day_drought  = day_end_drought - day_start_drought + 1
    day_all      = day_end - day_start_all + 1
    case_sum     = len(fcables)
    fw           = pd.DataFrame(np.zeros((day_drought+day_all)*case_sum),columns=['fwsoil'])
    fw['year']   = [''] * ((day_drought+day_all)*case_sum)
    fw['exp']    = [''] * ((day_drought+day_all)*case_sum)

    s = 0

    for case_num in np.arange(case_sum):

        cable = nc.Dataset(fcables[case_num], 'r')
        Time  = nc.num2date(cable.variables['time'][:],cable.variables['time'].units)

        Fwsoil          = pd.DataFrame(cable.variables['Fwsoil'][:,0,0],columns=['fwsoil'])
        Fwsoil['dates'] = Time
        Fwsoil          = Fwsoil.set_index('dates')
        Fwsoil          = Fwsoil.resample("D").agg('mean')
        Fwsoil.index    = Fwsoil.index - pd.datetime(2011,12,31)
        Fwsoil.index    = Fwsoil.index.days

        e  = s+day_drought

        fw['fwsoil'].iloc[s:e] = Fwsoil[np.all([Fwsoil.index >= day_start_drought,
                                 Fwsoil.index <=day_end_drought],axis=0)]['fwsoil'].values
        fw['year'].iloc[s:e]   = ['drought'] * day_drought
        fw['exp'].iloc[s:e]    = [ case_labels[case_num]] * day_drought
        s  = e
        e  = s+day_all
        fw['fwsoil'].iloc[s:e] = Fwsoil[np.all([Fwsoil.index >= day_start_all,
                                 Fwsoil.index <=day_end],axis=0)]['fwsoil'].values
        fw['year'].iloc[s:e]   = ['all'] * day_all
        fw['exp'].iloc[s:e]    = [ case_labels[case_num]] * day_all
        s  =  e

    # seaborn
    #sns.color_palette("Set2", 8)
    #flatui = ["orange", "dodgerblue"]
    #aaa = sns.set_palette(flatui)
    sns.boxplot(x="exp", y="fwsoil", hue="year", data=fw, palette="BrBG",
                order=case_labels,  width=0.7, hue_order=['drought','all'],
                ax=ax1, showfliers=False, color=almost_black) # palette="Set2",

    ax1.set_ylabel("$β$")
    ax1.set_xlabel("")
    ax1.axis('tight')
    ax1.set_ylim(0.,1.1)
    ax1.axhline(y=np.mean(fw[np.all([fw.year=='drought',fw.exp=='Ctl'],axis=0)]['fwsoil'].values),
                c=almost_black, ls="--")
    ax1.legend(loc='best', frameon=False)
    ax1.text(0.02, 0.95, '(a)', transform=ax1.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    print("***********************")
    # case_labels = ["Ctl", "Sres", "Watr", "Hi-Res-1", "Hi-Res-2", "Opt",  "β-hvrd",  "β-exp" ]
    print("median of Ctl is %f" % np.median(fw[np.all([fw.year=='all',fw.exp=='Ctl'],axis=0)]['fwsoil'].values))
    print("median of Sres is %f" % np.median(fw[np.all([fw.year=='all',fw.exp=='Sres'],axis=0)]['fwsoil'].values))
    print("median of Watr is %f" % np.median(fw[np.all([fw.year=='all',fw.exp=='Watr'],axis=0)]['fwsoil'].values))
    print("median of Hi-Res-1 is %f" % np.median(fw[np.all([fw.year=='all',fw.exp=='Hi-Res-1'],axis=0)]['fwsoil'].values))
    print("median of Hi-Res-2 is %f" % np.median(fw[np.all([fw.year=='all',fw.exp=='Hi-Res-2'],axis=0)]['fwsoil'].values))
    print("median of Opt is %f" % np.median(fw[np.all([fw.year=='all',fw.exp=='Opt'],axis=0)]['fwsoil'].values))
    print("median of β-hvrd is %f" % np.median(fw[np.all([fw.year=='all',fw.exp=='β-hvrd'],axis=0)]['fwsoil'].values))
    print("median of β-exp is %f" % np.median(fw[np.all([fw.year=='all',fw.exp=='β-exp'],axis=0)]['fwsoil'].values))
    print("***********************")

    #colors = cm.Set3(np.arange(0,len(case_labels)))
    colors = cm.tab20(np.arange(0,len(case_labels)))
    # ============================= boxplot ===================================
    for case_num in np.arange(len(fcables)):
        SM  = read_cable_SM(fcables[case_num], layers[case_num])
        fw  = read_cable_var(fcables[case_num], "Fwsoil")
        # print(SM)

        # theta_1.5m : using root zone soil moisture
        if layers[case_num] == "6":
            sm =(  SM.iloc[:,0]*0.022 + SM.iloc[:,1]*0.058 \
                 + SM.iloc[:,2]*0.154 + SM.iloc[:,3]*0.409 \
                 + SM.iloc[:,4]*(1.5-0.022-0.058-0.154-0.409) )/1.5
        elif layers[case_num] == "31uni":
            sm = SM.iloc[:,0:10].mean(axis = 1)

        # theta_all : using whole soil column soil moisture
        # if layers[case_num] == "6":
        #     sm =(  SM.iloc[:,0]*0.022 + SM.iloc[:,1]*0.058 \
        #          + SM.iloc[:,2]*0.154 + SM.iloc[:,3]*0.409 \
        #          + SM.iloc[:,4]*1.085 + SM.iloc[:,5]*2.872 )/4.6
        # elif layers[case_num] == "31uni":
        #     sm = SM.iloc[:,:].mean(axis = 1)

        ax2.scatter(sm, fw,  s=1., marker='o', alpha=0.45, c=colors[case_num],label=case_labels[case_num])

    ax2.set_xlim(0.08,0.405)
    ax2.set_ylim(0.0,1.05)
    ax2.set_ylabel("$β$")
    ax2.set_xlabel("$θ$$_{1.5m}$ (m$^{3}$ m$^{-3}$)")
    #ax2.set_xlabel("$θ$ (m$^{3}$ m$^{-3}$)")
    #ax2.legend(numpoints=1, loc='lower right', frameon=False)
    ax2.legend(numpoints=1, loc='best', frameon=False)
    ax2.text(0.02, 0.95, '(b)', transform=ax2.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    #plt.setp(ax2.get_yticklabels(), visible=False)

    fig.savefig("./plots/EucFACE_Fwsoil_boxplot_SM" , bbox_inches='tight', pad_inches=0.1)
