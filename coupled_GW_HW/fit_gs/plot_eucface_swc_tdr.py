#!/usr/bin/env python

"""
Calculate tdr SM, fwsoil, and fluxes

Include functions :

    plot_tdr
    plot_Fwsoil
    plot_ET
    plot_GPP
    plot_Rain
    plot_Rain_Fwsoil
    plot_ET_3

"""
__author__ = "MU Mengyuan"
__version__ = "2020-03-10"

import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import datetime as dt
import netCDF4 as nc
from matplotlib import cm
from matplotlib import ticker
from scipy.interpolate import griddata
import scipy.stats as stats
from sklearn.metrics import mean_squared_error
from plot_eucface_get_var import *

def plot_tdr(fcable, case_name, ring, layer):

    subset = read_obs_swc_tdr(ring)

# _________________________ CABLE ___________________________
    cable = nc.Dataset(fcable, 'r')
    Time  = nc.num2date(cable.variables['time'][:],cable.variables['time'].units)
    SoilMoist = pd.DataFrame(cable.variables['SoilMoist'][:,0,0,0], columns=['SoilMoist'])

    if layer == "6":
        SoilMoist['SoilMoist'] = ( cable.variables['SoilMoist'][:,0,0,0]*0.022 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058 \
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154 \
                                 + cable.variables['SoilMoist'][:,3,0,0]*(0.25-0.022-0.058-0.154) )/0.25
    elif layer == "31uni":
        SoilMoist['SoilMoist'] = ( cable.variables['SoilMoist'][:,0,0,0]*0.15 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.10 )/0.25

    SoilMoist['dates'] = Time
    SoilMoist = SoilMoist.set_index('dates')
    SoilMoist = SoilMoist.resample("D").agg('mean')
    SoilMoist.index = SoilMoist.index - pd.datetime(2011,12,31)
    SoilMoist.index = SoilMoist.index.days
    SoilMoist = SoilMoist.sort_values(by=['dates'])

    swilt = np.zeros(len(SoilMoist))
    sfc = np.zeros(len(SoilMoist))
    ssat = np.zeros(len(SoilMoist))

    if layer == "6":
        swilt[:] = ( cable.variables['swilt'][0]*0.022 + cable.variables['swilt'][1]*0.058 \
                   + cable.variables['swilt'][2]*0.154 + cable.variables['swilt'][3]*(0.25-0.022-0.058-0.154) )/0.25
        sfc[:] = ( cable.variables['sfc'][0]*0.022   + cable.variables['sfc'][1]*0.058 \
                   + cable.variables['sfc'][2]*0.154 + cable.variables['sfc'][3]*(0.25-0.022-0.058-0.154) )/0.25
        ssat[:] = ( cable.variables['ssat'][0]*0.022 + cable.variables['ssat'][1]*0.058 \
                   + cable.variables['ssat'][2]*0.154+ cable.variables['ssat'][3]*(0.25-0.022-0.058-0.154) )/0.25
    elif layer == "31uni":
        swilt[:] =(cable.variables['swilt'][0]*0.15 + cable.variables['swilt'][1]*0.10 )/0.25
        sfc[:] =(cable.variables['sfc'][0]*0.15 + cable.variables['sfc'][1]*0.10 )/0.25
        ssat[:] =(cable.variables['ssat'][0]*0.15 + cable.variables['ssat'][1]*0.10 )/0.25

# ____________________ Plot obs _______________________
    fig = plt.figure(figsize=[9,4])

    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color'] = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor'] = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    ax = fig.add_subplot(111)

    x   = SoilMoist.index

    ax.plot(subset.index, subset.values,   c="green", lw=1.0, ls="-", label="tdr")
    ax.plot(x, SoilMoist.values,c="orange", lw=1.0, ls="-", label="swc")
    '''
    tmp1 = SoilMoist['SoilMoist'].loc[SoilMoist.index.isin(subset.index)]
    tmp2 = subset.loc[subset.index.isin(SoilMoist.index)]
    mask = np.isnan(tmp2)
    print(mask)
    tmp1 = tmp1[mask == False]
    tmp2 = tmp2[mask == False]

    cor_tdr = stats.pearsonr(tmp1,tmp2)
    mse_tdr = mean_squared_error(tmp2, tmp1)
    ax.set_title("r = % 5.3f , MSE = % 5.3f" %(cor_tdr[0], np.sqrt(mse_tdr)))
    print("-----------------------------------------------")
    print(mse_tdr)
    '''
    ax.plot(x, swilt,           c="black", lw=1.0, ls="-", label="swilt")
    ax.plot(x, sfc,             c="black", lw=1.0, ls="-.", label="sfc")
    ax.plot(x, ssat,            c="black", lw=1.0, ls=":", label="ssat")

    cleaner_dates = ["2013","2014","2015","2016","2017","2018","2019"]
    xtickslocs    = [367,732,1097,1462,1828,2193,2558]

    plt.setp(ax.get_xticklabels(), visible=True)
    ax.set(xticks=xtickslocs, xticklabels=cleaner_dates) ####
    ax.set_ylabel("VWC (m3/m3)")
    ax.axis('tight')
    ax.set_ylim(0,0.5)
    ax.set_xlim(367,2739)
    ax.legend()

    fig.savefig("../plots/EucFACE_tdr_%s_%s" % (os.path.basename(case_name).split("/")[-1], ring), bbox_inches='tight', pad_inches=0.1)

def plot_Fwsoil(fcbl_def, fcbl_fw_def, fcbl_fw_hie, ring):

    fw1 = read_cable_var(fcbl_def, "Fwsoil")
    fw2 = read_cable_var(fcbl_fw_def, "Fwsoil")
    fw3 = read_cable_var(fcbl_fw_hie, "Fwsoil")

    fig = plt.figure(figsize=[15,10])

    ax  = fig.add_subplot(111)

    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color'] = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor'] = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    x = fw1.index

    ax.plot(x, fw1["cable"],   c="orange", lw=1.0, ls="-", label="Default_fw-std")
    ax.plot(x, fw2["cable"],   c="blue", lw=1.0, ls="-", label="Best_fw-std")
    ax.plot(x, fw3["cable"],   c="forestgreen", lw=1.0, ls="-", label="Best_fw-hie")

    cleaner_dates = ["2013","2014","2015","2016","2017","2018","2019"]
    xtickslocs    = [367,732,1097,1462,1828,2193,2558]

    plt.setp(ax.get_xticklabels(), visible=True)
    ax.set(xticks=xtickslocs, xticklabels=cleaner_dates)
    ax.set_ylabel("β")
    ax.axis('tight')
    ax.set_ylim(0.,1.1)
    ax.set_xlim(367,2739)
    ax.legend()

    fig.savefig("../plots/EucFACE_fwsoil_comp_%s" % ring, bbox_inches='tight', pad_inches=0.1)

def plot_check_ET(fcables, case_labels, ring):

    '''
    Check ET
    '''

    # ======================= PLOTTING  ==========================
    fig = plt.figure(figsize=[9,7])

    fig.subplots_adjust(hspace=0.15)
    fig.subplots_adjust(wspace=0.05)

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
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']  = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor']  = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')

    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    case_sum = len(fcables)
    colors   = cm.tab20(np.linspace(0,1,case_sum))

    # set x-axis values
    cleaner_dates1 = ["2013","2014","2015","2016","2017","2018","2019"]
    xtickslocs1    = [   367,   732,  1097,  1462,  1828,  2193, 2558 ]

    # ========================= ET FLUX  ============================
    # ===== Obs   =====
    subs_Esoil = read_obs_esoil(ring)
    subs_Trans = read_obs_trans(ring)

    ax1.scatter(subs_Trans.index, subs_Trans['obs'].rolling(window=3).mean(), marker='o', c='',edgecolors='red', s = 2., label="Obs")
    ax2.scatter(subs_Esoil.index, subs_Esoil['obs'].rolling(window=3).mean(), marker='o', c='',edgecolors='red', s = 2., label="Obs")
    # .rolling(window=3).mean()

    # ===== CABLE =====
    for i in np.arange(case_sum):
        cable = nc.Dataset(fcables[i], 'r')
        Time  = nc.num2date(cable.variables['time'][:],cable.variables['time'].units)

        TVeg = pd.DataFrame(cable.variables['TVeg'][:,0,0],columns=['TVeg'])
        TVeg = TVeg*1800.
        TVeg['dates'] = Time
        TVeg = TVeg.set_index('dates')
        TVeg = TVeg.resample("D").agg('sum')
        TVeg.index = TVeg.index - pd.datetime(2011,12,31)
        TVeg.index = TVeg.index.days

        ESoil = pd.DataFrame(cable.variables['ESoil'][:,0,0],columns=['ESoil'])
        ESoil = ESoil*1800.
        ESoil['dates'] = Time
        ESoil = ESoil.set_index('dates')
        ESoil = ESoil.resample("D").agg('sum')
        ESoil.index = ESoil.index - pd.datetime(2011,12,31)
        ESoil.index = ESoil.index.days

        x = TVeg.index

        #ax1.plot(x, TVeg['TVeg'].rolling(window=3).mean(),   c=colors[i], lw=1.0, ls="-", label=case_labels[i]) #
        #ax2.plot(x, ESoil['ESoil'].rolling(window=3).mean(), c=colors[i], lw=1.0, ls="-", label=case_labels[i]) #.rolling(window=7).mean()

        print(len(subs_Trans))
        print(len(TVeg[TVeg.index.isin(subs_Trans.index)]['TVeg']))
        print(len(x))
        print(len(TVeg['TVeg']))

        ax1.scatter(subs_Trans.index, TVeg[TVeg.index.isin(subs_Trans.index)]['TVeg'].rolling(window=3).mean(),
                    marker='o', c='', edgecolors=colors[i], s = 2., label=case_labels[i])

        ax2.scatter(subs_Esoil.index, ESoil[ESoil.index.isin(subs_Esoil.index)]['ESoil'].rolling(window=3).mean(),
                    marker='o', c='', edgecolors=colors[i], s = 2., label=case_labels[i])

        subs_Esoil['obs'].fillna(0., inplace=True)

        print(ESoil[np.all([ESoil.index.isin(subs_Esoil.index),np.isnan(subs_Esoil.values) == False],axis=0)]['ESoil'])
        print(subs_Esoil[np.isnan(subs_Esoil.values) == False]['obs'])
        print(stats.pearsonr(subs_Esoil[np.isnan(subs_Esoil.values) == False]['obs'],
              ESoil[np.all([ESoil.index.isin(subs_Esoil.index),np.isnan(subs_Esoil.values) == False],axis=0)]['ESoil']))

    ax1.text(0.02, 0.95, '(a)', transform=ax1.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax2.text(0.02, 0.95, '(b)', transform=ax2.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    # this order of the setting can affect plot x & y axis
    plt.setp(ax1.get_xticklabels(), visible=True)
    ax1.set(xticks=xtickslocs1, xticklabels=cleaner_dates1) ####
    ax1.set_ylabel("T (mm d$^{-1}$)")
    ax1.axis('tight')
    ax1.set_ylim(0.,3.5)
    ax1.set_xlim(367,1097)
    ax1.legend(loc='best', frameon=False)

    # this order of the setting can affect plot x & y axis
    plt.setp(ax2.get_xticklabels(), visible=True)
    ax2.set(xticks=xtickslocs1, xticklabels=cleaner_dates1) ####
    ax2.set_ylabel("Es (mm d$^{-1}$)")
    ax2.axis('tight')
    ax2.set_ylim(0.,4.5)
    ax2.set_xlim(367,1097)
    ax2.legend(loc='best', frameon=False)

    fig.savefig("../plots/EucFACE_Check_ET_%s.png" % ring, bbox_inches='tight', pad_inches=0.1)

def plot_ET(fcables, ring, case_labels):

    '''
    plot Trans and ESoil during 2014-1-15 ~ 2014-1-20 heatwave
    '''

    fig = plt.figure(figsize=[9,5])
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color'] = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor'] = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    ax  = fig.add_subplot(111)

    subs_Esoil = read_obs_esoil(ring)
    subs_Trans = read_obs_trans(ring)
    case_sum = len(fcables)

    colors = cm.tab20(np.linspace(0,1,case_sum))
    ax.scatter(subs_Trans.index, subs_Trans['obs'], marker='o', c='',edgecolors="green", s = 4., label="Trans Obs") # subs['EfloorPred'] 'blue'
    ax.scatter(subs_Esoil.index, subs_Esoil['obs'], marker='o', c='',edgecolors="red", s = 4., label="ESoil Obs") # subs['EfloorPred'] 'red'
    for case_num in np.arange(case_sum):
        TVeg  = read_cable_var(fcables[case_num], "TVeg")
        ESoil = read_cable_var(fcables[case_num], "ESoil")
        x = TVeg.index

        ax.plot(x, TVeg['cable'],     c=colors[case_num], lw=1.0, ls="-", label=case_labels[case_num]) #.rolling(window=5).mean() .rolling(window=7).mean()
        ax.plot(x, ESoil['cable'],    c=colors[case_num], lw=1.0, ls="-.") #.rolling(window=7).mean()

    #cleaner_dates = ["2013","2014","2015","2016","2017","2018","2019"] MMY
    #xtickslocs    = [367,732,1097,1462,1828,2193,2558] MMY
    cleaner_dates = ["2014-1-13","1-14","1-15","1-16","1-17","1-18","1-19","1-20"]
    xtickslocs    = [744,745,746,747,748,749,750,751]

    plt.setp(ax.get_xticklabels(), visible=True)
    ax.set(xticks=xtickslocs, xticklabels=cleaner_dates) ####
    ax.set_ylabel("Trans, Esoil (mm d$^{-1}$)")
    ax.axis('tight')
    ax.set_ylim(0.,3.0)
    ax.set_xlim(744,751)#(367,2739) # MMY
    # ax.set_xlim(367,1098) MMY
    ax.legend()

    fig.savefig("../plots/EucFACE_ET_%s" % ring, bbox_inches='tight', pad_inches=0.1)

def plot_GPP(fcables, ring, case_labels):

    '''
    plot GPP
    '''

    fig = plt.figure(figsize=[9,5])
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color'] = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor'] = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    ax  = fig.add_subplot(111)

    case_sum = len(fcables)

    colors = cm.tab20(np.linspace(0,1,case_sum))

    umol_2_gC = 12.0107 * 1.0E-6
    for case_num in np.arange(case_sum):
        GPP  = read_cable_GPP_year(fcables[case_num], "GPP")

        x = GPP.index

        ax.plot(x, GPP['cable']*umol_2_gC,     c=colors[case_num], lw=1.0, ls="-", label=case_labels[case_num]) #.rolling(window=5).mean() .rolling(window=7).mean()
        # .rolling(window=30).mean()*umol_2_gC

    #cleaner_dates = ["2013","2014","2015","2016","2017","2018","2019"]
    #xtickslocs    = [367,732,1097,1462,1828,2193,2558]

    #plt.setp(ax.get_xticklabels(), visible=True)
    #ax.set(xticks=xtickslocs, xticklabels=cleaner_dates) ####
    #ax.set_ylabel("GPP (umol m$^{-2}$)")
    ax.set_ylabel("GPP (gC m$^{-2}$ yr$^{-1}$)")
    ax.axis('tight')
    #ax.set_ylim(0.,3.0)
    #ax.set_xlim(367,2739) # MMY
    # ax.set_xlim(367,1098) MMY
    ax.legend()


    fig.savefig("../plots/EucFACE_GPP_%s" % ring, bbox_inches='tight', pad_inches=0.1)

def plot_Rain(fcable, case_name, ring):

    Rain  = read_cable_var(fcable, "Rainf")
    fig   = plt.figure(figsize=[15,10])

    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color'] = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor'] = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    ax    = fig.add_subplot(111)
    x     = Rain.index
    width = 1.

    ax.bar(x, Rain['cable'], width, color='royalblue', label='Obs')

    cleaner_dates = ["2013","2014","2015","2016","2017","2018","2019"]
    xtickslocs    = [367,732,1097,1462,1828,2193,2558]

    plt.setp(ax.get_xticklabels(), visible=True)
    ax.set(xticks=xtickslocs, xticklabels=cleaner_dates) ####
    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")
    ax.set_ylabel("Rain (mm/day)")
    ax.axis('tight')
    ax.set_ylim(0.,150.)
    ax.set_xlim(367,2739)

    fig.savefig("../plots/EucFACE_Rainfall", bbox_inches='tight', pad_inches=0.1)

def plot_Rain_Fwsoil(fcbl_def, fcbl_fw_def, fcbl_fw_hie, ring):

    fw1 = read_cable_var(fcbl_def, "Fwsoil")
    fw2 = read_cable_var(fcbl_fw_def, "Fwsoil")
    fw3 = read_cable_var(fcbl_fw_hie, "Fwsoil")
    Rain= read_cable_var(fcbl_def, "Rainf")

    fig = plt.figure(figsize=[15,10])

    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color'] = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor'] = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    ax1  = fig.add_subplot(211)
    ax2  = fig.add_subplot(212)

    x    = Rain.index
    width= 1.

    ax1.bar(x, Rain['cable'], width, color='royalblue', label='Obs')

    ax2.plot(x, fw1['cable'],   c="orange", lw=1.0, ls="-", label="Default_fw-std")
    ax2.plot(x, fw2['cable'],   c="blue", lw=1.0, ls="-", label="Best_fw-std")
    ax2.plot(x, fw3['cable'],   c="forestgreen", lw=1.0, ls="-", label="Best_fw-hie")

    cleaner_dates = ["2013","2014","2015","2016","2017","2018","2019"]
    xtickslocs    = [367,732,1097,1462,1828,2193,2558]

    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set(xticks=xtickslocs, xticklabels=cleaner_dates) ####
    ax1.yaxis.tick_left()
    ax1.yaxis.set_label_position("left")
    ax1.set_ylabel("Rain (mm/day)")
    ax1.axis('tight')
    ax1.set_ylim(0.,150.)
    ax1.set_xlim(367,2739)#,1098)

    plt.setp(ax2.get_xticklabels(), visible=True)
    ax2.set(xticks=xtickslocs, xticklabels=cleaner_dates)
    ax2.set_ylabel("β")
    ax2.axis('tight')
    ax2.set_ylim(0.,1.1)
    ax2.set_xlim(367,2739)#,1098)
    ax2.legend()

    fig.savefig("../plots/EucFACE_Rain_Fwsoil_%s" % ring, bbox_inches='tight', pad_inches=0.1)

def plot_ET_3(fctl, flit, fbest, ring):

    subs_Esoil = read_obs_esoil(ring)
    subs_Trans = read_obs_trans(ring)

    TVeg_ctl   = read_cable_var(fctl, "TVeg")
    ESoil_ctl  = read_cable_var(fctl, "ESoil")

    TVeg_lit   = read_cable_var(flit, "TVeg")
    ESoil_lit  = read_cable_var(flit, "ESoil")

    TVeg_best  = read_cable_var(fbest, "TVeg")
    ESoil_best = read_cable_var(fbest, "ESoil")

    fig = plt.figure(figsize=[9,12])
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color'] = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor'] = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    ax1  = fig.add_subplot(311)
    ax2  = fig.add_subplot(312)
    ax3  = fig.add_subplot(313)

    x = TVeg_ctl.index

    ax1.plot(x, TVeg_ctl['cable'].rolling(window=7).mean(),     c="green", lw=1.0, ls="-", label="Trans") #.rolling(window=5).mean() .rolling(window=7).mean()
    ax1.plot(x, ESoil_ctl['cable'].rolling(window=7).mean(),    c="orange", lw=1.0, ls="-", label="ESoil") #.rolling(window=7).mean()
    ax1.scatter(subs_Trans.index, subs_Trans['obs'], marker='o', c='',edgecolors="green", s = 4., label="Trans Obs") # subs['EfloorPred'] 'blue'
    ax1.scatter(subs_Esoil.index, subs_Esoil['obs'], marker='o', c='',edgecolors="orange", s = 4., label="ESoil Obs") # subs['EfloorPred'] 'red'

    ax2.plot(x, TVeg_lit['cable'].rolling(window=7).mean(),     c="green", lw=1.0, ls="-", label="Trans") #.rolling(window=5).mean() .rolling(window=7).mean()
    ax2.plot(x, ESoil_lit['cable'].rolling(window=7).mean(),    c="orange", lw=1.0, ls="-", label="ESoil") #.rolling(window=7).mean()
    ax2.scatter(subs_Trans.index, subs_Trans['obs'], marker='o', c='',edgecolors="green", s = 4., label="Trans Obs") # subs['EfloorPred'] 'blue'
    ax2.scatter(subs_Esoil.index, subs_Esoil['obs'], marker='o', c='',edgecolors="orange", s = 4., label="ESoil Obs") # subs['EfloorPred'] 'red'

    ax3.plot(x, TVeg_best['cable'].rolling(window=7).mean(),     c="green", lw=1.0, ls="-", label="Trans") #.rolling(window=5).mean() .rolling(window=7).mean()
    ax3.plot(x, ESoil_best['cable'].rolling(window=7).mean(),    c="orange", lw=1.0, ls="-", label="ESoil") #.rolling(window=7).mean()
    ax3.scatter(subs_Trans.index, subs_Trans['obs'], marker='o', c='',edgecolors="green", s = 4., label="Trans Obs") # subs['EfloorPred'] 'blue'
    ax3.scatter(subs_Esoil.index, subs_Esoil['obs'], marker='o', c='',edgecolors="orange", s = 4., label="ESoil Obs") # subs['EfloorPred'] 'red'

    cleaner_dates = ["2013","2014","2015","2016","2017","2018","2019"]
    xtickslocs    = [367,732,1097,1462,1828,2193,2558]

    plt.setp(ax1.get_xticklabels(), visible=True)
    ax1.set(xticks=xtickslocs, xticklabels=cleaner_dates) ####
    ax1.set_ylabel("Trans, Esoil ($mm d^{-1}$)")
    ax1.axis('tight')
    ax1.set_ylim(0.,4.0)
    ax1.set_xlim(367,1098)
    ax1.legend()

    plt.setp(ax2.get_xticklabels(), visible=True)
    ax2.set(xticks=xtickslocs, xticklabels=cleaner_dates) ####
    ax2.set_ylabel("Trans, Esoil ($mm d^{-1}$)")
    ax2.axis('tight')
    ax2.set_ylim(0.,4.0)
    ax2.set_xlim(367,1098)

    plt.setp(ax3.get_xticklabels(), visible=True)
    ax3.set(xticks=xtickslocs, xticklabels=cleaner_dates) ####
    ax3.set_ylabel("Trans, Esoil ($mm d^{-1}$)")
    ax3.axis('tight')
    ax3.set_ylim(0.,4.0)
    ax3.set_xlim(367,1098)

    fig.savefig("../plots/EucFACE_ET_ctl-lit-best" , bbox_inches='tight', pad_inches=0.1)
