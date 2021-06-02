#!/usr/bin/env python

"""
Get variable

include functions:

    read_cable_var
    read_cable_SM_one_clmn
    read_cable_SM
    read_obs_esoil
    read_obs_trans
    read_obs_swc_tdr
    read_obs_swc_neo
    read_SM_top_mid_bot
    read_obs_neo_top_mid_bot
    read_ET_SM_top_mid_bot

"""

__author__ = "MU Mengyuan"
__version__ = "2020-03-04"

import os
import sys
import numpy as np
import pandas as pd
import datetime as dt
import netCDF4 as nc
from scipy.interpolate import griddata

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
        var = var.resample("D").agg('sum')
    else:
        print("is here")
        var = var.resample("D").agg('mean')
    var.index = var.index - pd.datetime(2011,12,31)
    var.index = var.index.days
    var = var.sort_values(by=['Date'])

    return var

def read_cable_GPP_year(fcable, var_name):

    """
    read a var from CABLE output
    """

    print("carry on read_cable_var")
    cable = nc.Dataset(fcable, 'r')
    Time  = nc.num2date(cable.variables['time'][:],cable.variables['time'].units)
    var = pd.DataFrame(cable.variables[var_name][:,0,0]*1800., columns=['cable'])
    var['Date'] = Time
    var = var.set_index('Date')
    var = var.resample("Y").agg('sum')

    #var.index = var.index - pd.datetime(2011,12,31)
    #var.index = var.index.days
    var = var.sort_values(by=['Date'])

    return var

def read_cable_SM_one_clmn(fcable, layer):
    """
    Note: the SM here is for plotting profile, it has been turned into one column
          by SoilMoist = SoilMoist.stack()
    """
    cable = nc.Dataset(fcable, 'r')

    Time = nc.num2date(cable.variables['time'][:],cable.variables['time'].units)
    if layer == "6":
        SoilMoist = pd.DataFrame(cable.variables['SoilMoist'][:,:,0,0], columns=[1.1, 5.1, 15.7, 43.85, 118.55, 316.4])
    elif layer == "31uni":
        SoilMoist = pd.DataFrame(cable.variables['SoilMoist'][:,:,0,0], columns = \
                   [7.5,   22.5 , 37.5 , 52.5 , 67.5 , 82.5 , 97.5 , \
                    112.5, 127.5, 142.5, 157.5, 172.5, 187.5, 202.5, \
                    217.5, 232.5, 247.5, 262.5, 277.5, 292.5, 307.5, \
                    322.5, 337.5, 352.5, 367.5, 382.5, 397.5, 412.5, \
                    427.5, 442.5, 457.5 ])
    elif layer == "31exp":
        SoilMoist = pd.DataFrame(cable.variables['SoilMoist'][:,:,0,0], columns = \
                    [ 1.021985, 2.131912, 2.417723, 2.967358, 3.868759, 5.209868,\
                    7.078627, 9.562978, 12.75086, 16.73022, 21.58899, 27.41512,\
                    34.29655, 42.32122, 51.57708, 62.15205, 74.1341 , 87.61115,\
                    102.6711, 119.402 , 137.8918, 158.2283, 180.4995, 204.7933,\
                    231.1978, 259.8008, 290.6903, 323.9542, 359.6805, 397.9571,\
                    438.8719 ])
    elif layer == "31para":
        SoilMoist = pd.DataFrame(cable.variables['SoilMoist'][:,:,0,0], columns = \
                   [ 1.000014,  3.47101, 7.782496, 14.73158, 24.11537, 35.73098, \
                     49.37551, 64.84607, 81.93976, 100.4537, 120.185 , 140.9308, \
                     162.4881, 184.6541, 207.2259, 230.    , 252.7742, 275.346 , \
                     297.512 , 319.0693, 339.8151, 359.5464, 378.0603, 395.154 , \
                     410.6246, 424.2691, 435.8847, 445.2685, 452.2176, 456.5291, \
                     459.0001 ])

    SoilMoist['dates'] = Time
    SoilMoist = SoilMoist.set_index('dates')
    SoilMoist = SoilMoist.resample("D").agg('mean')
    SoilMoist.index = SoilMoist.index - pd.datetime(2011,12,31)
    SoilMoist.index = SoilMoist.index.days
    SoilMoist = SoilMoist.stack() # turn multi-columns into one-column
    SoilMoist = SoilMoist.reset_index() # remove index 'dates'
    SoilMoist = SoilMoist.rename(index=str, columns={"level_1": "Depth"})
    SoilMoist = SoilMoist.sort_values(by=['Depth','dates'])
    # rename columns level_1 to Depth
    #SoilMoist = SoilMoist.set_index('Depth')
    print("----------------------------")
    print(SoilMoist)

    return SoilMoist

def read_cable_SM(fcable, layer):
    """
    Note: the SM here is a multi column dataframe, it doesn't aim to plot profile plot
    """
    cable = nc.Dataset(fcable, 'r')

    Time = nc.num2date(cable.variables['time'][:],cable.variables['time'].units)
    if layer == "6":
        SoilMoist = pd.DataFrame(cable.variables['SoilMoist'][:,:,0,0], columns=[1.1, 5.1, 15.7, 43.85, 118.55, 316.4])
    elif layer == "31uni":
        SoilMoist = pd.DataFrame(cable.variables['SoilMoist'][:,:,0,0], columns = \
                   [7.5,   22.5 , 37.5 , 52.5 , 67.5 , 82.5 , 97.5 , \
                    112.5, 127.5, 142.5, 157.5, 172.5, 187.5, 202.5, \
                    217.5, 232.5, 247.5, 262.5, 277.5, 292.5, 307.5, \
                    322.5, 337.5, 352.5, 367.5, 382.5, 397.5, 412.5, \
                    427.5, 442.5, 457.5 ])
    elif layer == "31exp":
        SoilMoist = pd.DataFrame(cable.variables['SoilMoist'][:,:,0,0], columns = \
                    [ 1.021985, 2.131912, 2.417723, 2.967358, 3.868759, 5.209868,\
                    7.078627, 9.562978, 12.75086, 16.73022, 21.58899, 27.41512,\
                    34.29655, 42.32122, 51.57708, 62.15205, 74.1341 , 87.61115,\
                    102.6711, 119.402 , 137.8918, 158.2283, 180.4995, 204.7933,\
                    231.1978, 259.8008, 290.6903, 323.9542, 359.6805, 397.9571,\
                    438.8719 ])
    elif layer == "31para":
        SoilMoist = pd.DataFrame(cable.variables['SoilMoist'][:,:,0,0], columns = \
                   [ 1.000014,  3.47101, 7.782496, 14.73158, 24.11537, 35.73098, \
                     49.37551, 64.84607, 81.93976, 100.4537, 120.185 , 140.9308, \
                     162.4881, 184.6541, 207.2259, 230.    , 252.7742, 275.346 , \
                     297.512 , 319.0693, 339.8151, 359.5464, 378.0603, 395.154 , \
                     410.6246, 424.2691, 435.8847, 445.2685, 452.2176, 456.5291, \
                     459.0001 ])

    SoilMoist['dates'] = Time
    SoilMoist = SoilMoist.set_index('dates')
    SoilMoist = SoilMoist.resample("D").agg('mean')
    SoilMoist.index = SoilMoist.index - pd.datetime(2011,12,31)
    SoilMoist.index = SoilMoist.index.days

    return SoilMoist

def read_obs_esoil(ring):

    '''
    Using wuTP

    wuTP: understrorey evapotranspiration in mm/day estimated from changes
          in upper soil moisture
    EfloorPred: understrorey evapotranspiration in mm/day estimated from
          the nonlinear correlation with site potential evapotranspiration
    '''

    fobs_Esoil = "/srv/ccrc/data25/z5218916/data/Eucface_data/FACE_PACKAGE_HYDROMET_GIMENO_20120430-20141115/data/Gimeno_wb_EucFACE_underET.csv"

    est_esoil = pd.read_csv(fobs_Esoil, usecols = ['Ring','Date','wuTP'])
    est_esoil['Date'] = pd.to_datetime(est_esoil['Date'],format="%d/%m/%Y",infer_datetime_format=False)
    est_esoil['Date'] = est_esoil['Date'] - pd.datetime(2011,12,31)
    est_esoil['Date'] = est_esoil['Date'].dt.days
    est_esoil = est_esoil.sort_values(by=['Date'])
    # divide neo into groups
    if ring == 'amb':
       subs = est_esoil[(est_esoil['Ring'].isin(['R2','R3','R6'])) & (est_esoil.Date > 366)]
    elif ring == 'ele':
       subs = est_esoil[(est_esoil['Ring'].isin(['R1','R4','R5'])) & (est_esoil.Date > 366)]
    else:
       subs = est_esoil[(est_esoil['Ring'].isin([ring]))  & (est_esoil.Date > 366)]

    subs = subs.groupby(by=["Date"]).mean()
    subs['wuTP']   = subs['wuTP'].clip(lower=0.)
    subs['wuTP']   = subs['wuTP'].replace(0., float('nan'))
    subs = subs.rename({'wuTP' : 'obs'}, axis='columns')

    return subs

def read_obs_trans(ring):

    fobs_Trans = "/srv/ccrc/data25/z5218916/data/Eucface_data/FACE_PACKAGE_HYDROMET_GIMENO_20120430-20141115/data/Gimeno_wb_EucFACE_sapflow.csv"

    est_trans = pd.read_csv(fobs_Trans, usecols = ['Ring','Date','volRing'])
    est_trans['Date'] = pd.to_datetime(est_trans['Date'],format="%d/%m/%Y",infer_datetime_format=False)
    est_trans['Date'] = est_trans['Date'] - pd.datetime(2011,12,31)
    est_trans['Date'] = est_trans['Date'].dt.days
    est_trans = est_trans.sort_values(by=['Date'])
    # divide neo into groups
    if ring == 'amb':
       subs = est_trans[(est_trans['Ring'].isin(['R2','R3','R6'])) & (est_trans.Date > 366)]
    elif ring == 'ele':
       subs = est_trans[(est_trans['Ring'].isin(['R1','R4','R5'])) & (est_trans.Date > 366)]
    else:
       subs = est_trans[(est_trans['Ring'].isin([ring]))  & (est_trans.Date > 366)]

    subs = subs.groupby(by=["Date"]).mean()
    subs['volRing']   = subs['volRing'].clip(lower=0.)
    subs['volRing']   = subs['volRing'].replace(0., float('nan'))
    subs = subs.rename({'volRing' : 'obs'}, axis='columns')

    return subs

def read_obs_swc_tdr(ring):

    fobs   = "/srv/ccrc/data25/z5218916/data/Eucface_data/SM_2013-2019/eucSM1319_gap_filled.csv"
    tdr = pd.read_csv(fobs, usecols = ['Ring','Date','swc.tdr'])
    tdr['Date'] = pd.to_datetime(tdr['Date'],format="%d/%m/%Y",infer_datetime_format=False) # "%Y-%m-%d"
    tdr['Date'] = tdr['Date'] - pd.datetime(2011,12,31)
    tdr['Date'] = tdr['Date'].dt.days
    tdr = tdr.sort_values(by=['Date'])
    # divide neo into groups
    if ring == 'amb':
        subset = tdr[(tdr['Ring'].isin(['R2','R3','R6'])) & (tdr.Date > 366)]
    elif ring == 'ele':
        subset = tdr[(tdr['Ring'].isin(['R1','R4','R5'])) & (tdr.Date > 366)]
    else:
        subset = tdr[(tdr['Ring'].isin([ring]))  & (tdr.Date > 366)]

    subset = subset.groupby(by=["Date"]).mean()/100.
    subset = subset.rename({'swc.tdr' : 'obs'}, axis='columns')
    return subset

def read_obs_swc_neo(ring):

    fobs = "/srv/ccrc/data25/z5218916/data/Eucface_data/swc_at_depth/FACE_P0018_RA_NEUTRON_20120430-20190510_L1.csv"
    neo = pd.read_csv(fobs, usecols = ['Ring','Depth','Date','VWC'])
    # usecols : read specific columns from CSV

    # translate datetime
    neo['Date'] = pd.to_datetime(neo['Date'],format="%d/%m/%y",infer_datetime_format=False)
    #  unit='D', origin=pd.Timestamp('2012-01-01')

    # turn datetime64[ns] into timedelta64[ns] since 2011-12-31, e.g. 2012-1-1 as 1 days
    neo['Date'] = neo['Date'] - pd.datetime(2011,12,31)

    # extract days as integers from a timedelta64[ns] object
    neo['Date'] = neo['Date'].dt.days

    # sort by 'Date','Depth'
    neo = neo.sort_values(by=['Date','Depth'])

    print(neo['Depth'].unique())

    # divide neo into groups
    if ring == 'amb':
        subset = neo[neo['Ring'].isin(['R2','R3','R6'])]
    elif ring == 'ele':
        subset = neo[neo['Ring'].isin(['R1','R4','R5'])]
    else:
        subset = neo[neo['Ring'].isin([ring])]

    # calculate the mean of every group ( and unstack #.unstack(level=0)
    subset = subset.groupby(by=["Depth","Date"]).mean()

    # remove 'VWC'
    subset = subset.xs('VWC', axis=1, drop_level=True)

    # 'VWC' : key on which to get cross section
    # axis=1 : get cross section of column
    # drop_level=True : returns cross section without the multilevel index

    #neo_mean = np.transpose(neo_mean)

    return subset

def read_SM_top_mid_bot(fcable, ring, layer):
    """
    Read CABLE ET and oil moisture for top mid bot blocks used in metrics calculation

    """

    print(layer)
    cable = nc.Dataset(fcable, 'r')
    Time  = nc.num2date(cable.variables['time'][:],cable.variables['time'].units)
    cable_data = pd.DataFrame(cable.variables['TVeg'][:,0,0]*1800., columns=['TVeg'])

    if layer == "6":
        cable_data['SM_top']  = (  cable.variables['SoilMoist'][:,0,0,0]*0.022 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058\
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154 \
                                 + cable.variables['SoilMoist'][:,3,0,0]*(0.3-0.022-0.058-0.154) )/0.3
        cable_data['SM_mid']  = (  cable.variables['SoilMoist'][:,3,0,0]*0.343 \
                                 + cable.variables['SoilMoist'][:,4,0,0]*(1.2-0.343) )/1.2
        cable_data['SM_bot']  = (  cable.variables['SoilMoist'][:,4,0,0]*(1.085-(1.2-0.343)) \
                                 + cable.variables['SoilMoist'][:,5,0,0]*2.872)/(4.6-1.5)
        cable_data['SM_all'] = (   cable.variables['SoilMoist'][:,0,0,0]*0.022 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058 \
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154 \
                                 + cable.variables['SoilMoist'][:,3,0,0]*0.409 \
                                 + cable.variables['SoilMoist'][:,4,0,0]*1.085 \
                                 + cable.variables['SoilMoist'][:,5,0,0]*2.872  )/4.6
        cable_data['SM_15m']  = (  cable.variables['SoilMoist'][:,0,0,0]*0.022 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058 \
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154 \
                                 + cable.variables['SoilMoist'][:,3,0,0]*0.409 \
                                 + cable.variables['SoilMoist'][:,4,0,0]*      \
                                 (1.5-0.022-0.058-0.154-0.409) )/1.5
    elif layer == "31uni":
        print("come in")
        cable_data['SM_top']  = (cable.variables['SoilMoist'][:,0,0,0]*0.15 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.15)/0.3
        cable_data['SM_mid']  = cable.variables['SoilMoist'][:,2,0,0]*0.15
        for i in np.arange(3,10):
            cable_data['SM_mid']  = cable_data['SM_mid'] + cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_mid']  = cable_data['SM_mid']/(1.5-0.3)

        cable_data['SM_bot']  = cable.variables['SoilMoist'][:,10,0,0]*0.15
        for i in np.arange(11,30):
            cable_data['SM_bot']  = cable_data['SM_bot'] + cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_bot']  = (cable_data['SM_bot'] + cable.variables['SoilMoist'][:,30,0,0]*0.1)/(4.6-1.5)

        cable_data['SM_all']  = cable.variables['SoilMoist'][:,30,0,0]*0.1
        for i in np.arange(0,30):
            cable_data['SM_all']  = cable_data['SM_all'] + cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_all'] = cable_data['SM_all']/4.6

        cable_data['SM_15m']  = cable.variables['SoilMoist'][:,0,0,0]*0.15
        for i in np.arange(1,10):
            cable_data['SM_15m']  = cable_data['SM_15m']+ cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_15m']  = cable_data['SM_15m']/1.5

    cable_data['dates'] = Time
    cable_data = cable_data.set_index('dates')
    cable_data = cable_data.resample("D").agg('mean')
    cable_data.index = cable_data.index - pd.datetime(2011,12,31)
    cable_data.index = cable_data.index.days
    cable_data = cable_data.sort_values(by=['dates'])

    return cable_data

def read_SM_top_mid_bot_hourly(fcable, ring, layer):

    """
    Read CABLE ET and oil moisture for top mid bot blocks used in metrics calculation

    """

    cable = nc.Dataset(fcable, 'r')
    Time  = nc.num2date(cable.variables['time'][:],cable.variables['time'].units)
    cable_data = pd.DataFrame(cable.variables['TVeg'][:,0,0]*1800., columns=['TVeg'])

    if layer == "6":
        cable_data['SM_top']  = (  cable.variables['SoilMoist'][:,0,0,0]*0.022 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058\
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154 \
                                 + cable.variables['SoilMoist'][:,3,0,0]*(0.3-0.022-0.058-0.154) )/0.3
        cable_data['SM_mid']  = (  cable.variables['SoilMoist'][:,3,0,0]*0.343 \
                                 + cable.variables['SoilMoist'][:,4,0,0]*(1.2-0.343) )/1.2
        cable_data['SM_bot']  = (  cable.variables['SoilMoist'][:,4,0,0]*(1.085-(1.2-0.343)) \
                                 + cable.variables['SoilMoist'][:,5,0,0]*2.872)/(4.6-1.5)
        cable_data['SM_all'] = (   cable.variables['SoilMoist'][:,0,0,0]*0.022 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058 \
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154 \
                                 + cable.variables['SoilMoist'][:,3,0,0]*0.409 \
                                 + cable.variables['SoilMoist'][:,4,0,0]*1.085 \
                                 + cable.variables['SoilMoist'][:,5,0,0]*2.872  )/4.6
        cable_data['SM_15m']  = (  cable.variables['SoilMoist'][:,0,0,0]*0.022 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058 \
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154 \
                                 + cable.variables['SoilMoist'][:,3,0,0]*0.409 \
                                 + cable.variables['SoilMoist'][:,4,0,0]*      \
                                 (1.5-0.022-0.058-0.154-0.409) )/1.5
    elif layer == "31uni":
        cable_data['SM_top']  = (cable.variables['SoilMoist'][:,0,0,0]*0.15 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.15)/0.3
        cable_data['SM_mid']  = cable.variables['SoilMoist'][:,2,0,0]*0.15
        for i in np.arange(3,10):
            cable_data['SM_mid']  = cable_data['SM_mid'] + cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_mid']  = cable_data['SM_mid']/(1.5-0.3)

        cable_data['SM_bot']  = cable.variables['SoilMoist'][:,10,0,0]*0.15
        for i in np.arange(11,30):
            cable_data['SM_bot']  = cable_data['SM_bot'] + cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_bot']  = (cable_data['SM_bot'] + cable.variables['SoilMoist'][:,30,0,0]*0.1)/(4.6-1.5)

        cable_data['SM_all']  = cable.variables['SoilMoist'][:,30,0,0]*0.1
        for i in np.arange(0,30):
            cable_data['SM_all']  = cable_data['SM_all'] + cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_all'] = cable_data['SM_all']/4.6

        cable_data['SM_15m']  = cable.variables['SoilMoist'][:,0,0,0]*0.15
        for i in np.arange(1,10):
            cable_data['SM_15m']  = cable_data['SM_15m']+ cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_15m']  = cable_data['SM_15m']/1.5

    cable_data['dates'] = Time
    cable_data = cable_data.set_index('dates')

    return cable_data

def read_obs_neo_top_mid_bot(ring):

    """
    Read neo soil moisture for top mid and bot soil blocks used for metrics calculation
    """
    fobs_neo = "/srv/ccrc/data25/z5218916/data/Eucface_data/swc_at_depth/FACE_P0018_RA_NEUTRON_20120430-20190510_L1.csv"
    neo = pd.read_csv(fobs_neo, usecols = ['Ring','Depth','Date','VWC'])
    neo['Date'] = pd.to_datetime(neo['Date'],format="%d/%m/%y",infer_datetime_format=False)
    neo['Date'] = neo['Date'] - pd.datetime(2011,12,31)
    neo['Date'] = neo['Date'].dt.days
    neo = neo.sort_values(by=['Date','Depth'])

    if ring == 'amb':
        subset = neo[neo['Ring'].isin(['R2','R3','R6'])]
    elif ring == 'ele':
        subset = neo[neo['Ring'].isin(['R1','R4','R5'])]
    else:
        subset = neo[neo['Ring'].isin([ring])]

    subset = subset.groupby(by=["Depth","Date"]).mean()
    subset = subset.xs('VWC', axis=1, drop_level=True)
    x     = subset.index.get_level_values(1).values
    y     = subset.index.get_level_values(0).values
    value = subset.values

    X     = subset[(25)].index.values[20:]
    Y     = np.arange(0.5,460,1)

    grid_X, grid_Y = np.meshgrid(X,Y)

    grid_data = griddata((x, y) , value, (grid_X, grid_Y), method='nearest')

    neo_data = pd.DataFrame(subset[(25)].index.values[20:], columns=['dates'])
    neo_data["SM_15m"] = np.mean(grid_data[0:150,:],axis=0)/100.
    neo_data["SM_all"] = np.mean(grid_data[:,:],axis=0)/100.
    neo_data["SM_top"] = np.mean(grid_data[0:30,:],axis=0)/100.
    neo_data["SM_mid"]  = np.mean(grid_data[30:150,:],axis=0)/100.
    neo_data["SM_bot"] = np.mean(grid_data[150:460,:],axis=0)/100.
    neo_data["WA_all"] = np.sum(grid_data[:,:]/100.*10.,axis=0)
    neo_data = neo_data.set_index('dates')
    return neo_data

def read_ET_SM_top_mid_bot(fcable, ring, layer):

    """
    Read CABLE ET and oil moisture for top mid bot blocks used in metrics calculation
    """
    cable = nc.Dataset(fcable, 'r')
    Time  = nc.num2date(cable.variables['time'][:],cable.variables['time'].units)
    cable_data = pd.DataFrame(cable.variables['TVeg'][:,0,0]*1800., columns=['TVeg'])
    cable_data['ESoil'] = cable.variables['ESoil'][:,0,0]*1800.
    cable_data['Evap'] = cable.variables['Evap'][:,0,0]*1800.

    if layer == "6":
        cable_data['SM_25cm'] = (  cable.variables['SoilMoist'][:,0,0,0]*0.022
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154
                                 + cable.variables['SoilMoist'][:,3,0,0]*(0.25-0.022-0.058-0.154) )/0.25
        cable_data['SM_15m'] = (   cable.variables['SoilMoist'][:,0,0,0]*0.022
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154
                                 + cable.variables['SoilMoist'][:,3,0,0]*0.409
                                 + cable.variables['SoilMoist'][:,4,0,0]*
                                   (1.5-0.022-0.058-0.154-0.409))/1.5
        cable_data['SM_all'] = (  cable.variables['SoilMoist'][:,0,0,0]*0.022 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058 \
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154 \
                                 + cable.variables['SoilMoist'][:,3,0,0]*0.409 \
                                 + cable.variables['SoilMoist'][:,4,0,0]*1.085 \
                                 + cable.variables['SoilMoist'][:,5,0,0]*2.872  )/4.6
        cable_data['SM_top']  = (  cable.variables['SoilMoist'][:,0,0,0]*0.022 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058\
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154 \
                                 + cable.variables['SoilMoist'][:,3,0,0]*(0.3-0.022-0.058-0.154) )/0.3
        cable_data['SM_mid']  = (  cable.variables['SoilMoist'][:,3,0,0]*0.343 \
                                 + cable.variables['SoilMoist'][:,4,0,0]*(1.2-0.343) )/1.2
        cable_data['SM_bot']  = (  cable.variables['SoilMoist'][:,4,0,0]*(1.085-(1.2-0.343)) \
                                 + cable.variables['SoilMoist'][:,5,0,0]*2.872)/(4.6-1.5)
        cable_data['WA_all'] = (   cable.variables['SoilMoist'][:,0,0,0]*0.022 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.058 \
                                 + cable.variables['SoilMoist'][:,2,0,0]*0.154 \
                                 + cable.variables['SoilMoist'][:,3,0,0]*0.409 \
                                 + cable.variables['SoilMoist'][:,4,0,0]*1.085 \
                                 + cable.variables['SoilMoist'][:,5,0,0]*2.872  )*1000.

    elif layer == "31uni":
        cable_data['SM_25cm'] = ( cable.variables['SoilMoist'][:,0,0,0]*0.15 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.10 )/0.25
        cable_data['SM_15m']  = cable.variables['SoilMoist'][:,0,0,0]*0.15
        for i in np.arange(1,10):
            cable_data['SM_15m']  = cable_data['SM_15m'] + cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_15m']  = cable_data['SM_15m']/1.5

        cable_data['SM_all']  = cable.variables['SoilMoist'][:,30,0,0]*0.1
        for i in np.arange(0,30):
            cable_data['SM_all']  = cable_data['SM_all'] + cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_all']  = cable_data['SM_all']/4.6

        cable_data['SM_top']  = (cable.variables['SoilMoist'][:,0,0,0]*0.15 \
                                 + cable.variables['SoilMoist'][:,1,0,0]*0.15)/0.3
        cable_data['SM_mid']  = cable.variables['SoilMoist'][:,2,0,0]*0.15
        for i in np.arange(3,10):
            cable_data['SM_mid']  = cable_data['SM_mid'] + cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_mid']  = cable_data['SM_mid']/(1.5-0.3)

        cable_data['SM_bot']  = cable.variables['SoilMoist'][:,10,0,0]*0.15
        for i in np.arange(11,30):
            cable_data['SM_bot']  = cable_data['SM_bot'] + cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['SM_bot']  = (cable_data['SM_bot'] + cable.variables['SoilMoist'][:,30,0,0]*0.1)/(4.6-1.5)

        cable_data['WA_all']  = cable.variables['SoilMoist'][:,30,0,0]*0.1
        for i in np.arange(0,30):
            cable_data['WA_all']  = cable_data['WA_all'] + cable.variables['SoilMoist'][:,i,0,0]*0.15
        cable_data['WA_all'] = cable_data['WA_all']*1000.

    cable_data['dates'] = Time
    cable_data = cable_data.set_index('dates')
    cable_data = cable_data.resample("D").agg('mean')

    # from mm / half-hour to mm/day
    cable_data['TVeg']  = cable_data['TVeg']*48
    cable_data['ESoil'] = cable_data['ESoil']*48
    cable_data['Evap']  = cable_data['Evap']*48

    cable_data.index = cable_data.index - pd.datetime(2011,12,31)
    cable_data.index = cable_data.index.days
    cable_data = cable_data.sort_values(by=['dates'])

    return cable_data
