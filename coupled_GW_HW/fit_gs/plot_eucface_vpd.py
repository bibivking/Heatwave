#!/usr/bin/env python

"""
draw Trans_cable vs Trans_obs classified by EucFACE_Trans_VPD

Include functions :
    plot_vpd
    best_fit
    qair_to_vpd
    read_cable_var
    read_obs_trans

"""
__author__ = "MU Mengyuan"
__version__ = "2020-06-15"

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

def plot_vpd(fcable, ring):

    # Read in data
    press  = read_cable_var(fcable, "PSurf")[:8760] # surface air pressure
    qair   = read_cable_var(fcable, "Qair")[:8760]  # Surface specific humidity
    tair   = read_cable_var(fcable, "Tair")[:8760]  # Surface air temperature
    Tcable = read_cable_var(fcable, "TVeg")[:8760]
    Tobs   = read_obs_trans(ring)

    print(press)
    print(qair)
    print(Tobs)
    VPD   = qair_to_vpd(qair, tair, press)
    print(VPD)

    # VPD <= 1
    Tcable1 =  Tcable[VPD < 1.]
    Tobs1  =   Tobs[VPD < 1.]

    # 1 <= VPD < 1.5
    Tcable2 =  Tcable[np.all([VPD < 2,VPD >= 1.],axis=0)]
    Tobs2  =   Tobs[np.all([VPD < 2,VPD >= 1.],axis=0)]

    # 1.5 <= VPD < 2
    Tcable3 =  Tcable[np.all([VPD < 3.,VPD >= 2],axis=0)]
    Tobs3  =   Tobs[np.all([VPD < 3.,VPD >= 2],axis=0)]

    # 2 < VPD <= 2.5
    Tcable4 =  Tcable[VPD >= 3.]
    Tobs4  =   Tobs[VPD >= 3.]

    # fit solution
    a, b = best_fit(Tobs['obs'][:], Tcable['cable'][:])
    print(Tcable.values)

    x    = np.arange(0,0.45,0.01)
    yfit = [a + b * xi for xi in x]

    # _____________ Make plot _____________
    fig = plt.figure(figsize=(6,5))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams["legend.markerscale"] = 2.0

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
    colors = cm.Set2(np.arange(0,4))

    ax1 = fig.add_subplot(111)

    colors = cm.viridis(np.arange(0,4))

    ax1.plot(x, x, lw=1.5, ls="-", c = almost_black, alpha=1.)
    ax1.plot(x, yfit, lw=1.5, ls="--", c = almost_black, alpha=1.)
    # ax1.scatter( Tobs1, Tcable1, marker='o', c="darkslategray",edgecolors="darkslategray", s = 4., label="<1" , alpha=1)
    # ax1.scatter( Tobs2, Tcable2, marker='o', c="steelblue"         ,edgecolors="steelblue", s = 4., label="1-2" , alpha=1)
    ax1.scatter( Tobs3, Tcable3,marker='o', c="limegreen"        ,edgecolors="limegreen", s = 4., label="2-3" , alpha=1)
    #ax1.scatter(Tcable4, Tobs4, marker='o', c="springgreen"       ,edgecolors="springgreen", s = 6., label="2-2.5", alpha=1)
    ax1.scatter( Tobs4, Tcable4,marker='o', c="gold"        ,edgecolors="gold", s = 4., label=">3" , alpha=1)
    #ax1.scatter(Tcable6, Tobs6, marker='o', c="yellow"       ,edgecolors="yellow", s = 4., label=">3", alpha=0.5)


    ax1.axis('tight')
    ax1.set_xticks([0,0.1,0.2,0.3,0.4])
    ax1.set_xticklabels([0,0.1,0.2,0.3,0.4])
    ax1.set_yticks([0,0.1,0.2,0.3,0.4])
    ax1.set_yticklabels([0,0.1,0.2,0.3,0.4])
    ax1.set_xlim(0.,0.38)
    ax1.set_ylim(0.,0.38)
    ax1.set_xlabel('Measured $E_{tr}$ (mm hr$^{-1}$)')
    ax1.set_ylabel('Modelled $E_{tr}$ (mm hr$^{-1}$)')
    ax1.text(0.82, 0.25, '$D$ (kPa)', c = almost_black, transform=ax1.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    # 0.42
    ax1.legend(loc='lower right', frameon=False)

    fig.savefig("./plots/EucFACE_Trans_VPD_D-2+" , bbox_inches='tight', pad_inches=0.1)

def best_fit(X, Y):

    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2
    # print(numer)
    # print(denum)
    b = numer / denum
    a = ybar - b * xbar
    # print(b)
    # print(a)
    print('best fit line: y = %2f + %2f x' %(a, b))

    return a, b

def qair_to_vpd(qair, tair, press):
    '''
    calculate vpd
    '''
    DEG_2_KELVIN = 273.15
    PA_TO_KPA = 0.001
    PA_TO_HPA = 0.01

    # convert back to Pa
    press /= PA_TO_HPA
    tair  -= DEG_2_KELVIN

    # saturation vapor pressure
    es = 100.0 * 6.112 * np.exp((17.67 * tair) / (243.5 + tair))

    # vapor pressure
    ea = (qair * press) / (0.622 + (1.0 - 0.622) * qair)

    vpd = (es - ea) * PA_TO_KPA
    vpd = np.where(vpd < 0.05, 0.05, vpd)

    return vpd

def read_cable_var(fcable, var_name):

    """
    read a var from CABLE output
    """

    print("carry on read_cable_var")
    cable = nc.Dataset(fcable, 'r')
    Time  = nc.num2date(cable.variables['time'][:],cable.variables['time'].units)
    if var_name in ["TVeg", "ESoil", "Rainf", "GPP"]:
        var = pd.DataFrame(cable.variables[var_name][:,0,0]*1800., columns=['cable'])
    else:
        var = pd.DataFrame(cable.variables[var_name][:,0,0], columns=['cable'])
    var['Date'] = Time
    var = var.set_index('Date')
    if var_name in ["TVeg", "ESoil", "Rainf", "GPP"]:
        var = var.resample("H").agg('sum')
    else:
        print("is here")
        var = var.resample("H").agg('mean')
    var = var.sort_values(by=['Date'])

    return var

def read_obs_trans(ring):

    fobs_Trans = "/srv/ccrc/data25/z5218916/data/Eucface_data/Jim_vpd_data/eucface13-16.csv"

    est_trans = pd.read_csv(fobs_Trans, usecols = ['Ring','DateTime','sap'])
    est_trans['Date'] = pd.to_datetime(est_trans['DateTime'],format="%Y-%m-%d %H:%M:%S",infer_datetime_format=False)
    est_trans = est_trans.sort_values(by=['Date'])
    # divide neo into groups
    if ring == 'amb':
       subs = est_trans[(est_trans['Ring'].isin(['R2','R3','R6']))]
    elif ring == 'ele':
       subs = est_trans[(est_trans['Ring'].isin(['R1','R4','R5']))]
    else:
       subs = est_trans[(est_trans['Ring'].isin([ring]))]

    subs = subs.groupby(by=["Date"]).mean()
    # subs['sap']   = subs['sap'].clip(lower=0.)
    # subs['sap']   = subs['sap'].replace(0., float('nan'))
    subs = subs.rename({'sap' : 'obs'}, axis='columns')

    return subs

if __name__ == "__main__":


    ring    = "amb"
    fcable  = "/srv/ccrc/data25/z5218916/cable/EucFACE/EucFACE_run/outputs/met_LAI-08_vrt_swilt-watr-ssat_hyds10_31uni_teuc_sres_watr/EucFACE_amb_out.nc"

    plot_vpd(fcable, ring)
