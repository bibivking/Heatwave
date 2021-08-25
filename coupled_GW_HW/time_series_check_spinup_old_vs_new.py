#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from convert_units import get_land_var_scale
from common_utils import *

def read_LIS_CABLE(flis,var_name,var_dim,latidxs=None,lonidxs=None):

    data = Dataset(flis, mode='r')

    # Replace _FillValues with NaNs:
    if var_dim == 3:
        if latidxs == None and lonidxs == None:
            var_nans = data.variables[var_name][:,:,:]
        else:
            var_nans = data.variables[var_name][:,latidxs[0]:latidxs[1],lonidxs[0]:lonidxs[1]]
    elif var_dim == 4:
        if latidxs == None and lonidxs == None:
            var_nans = data.variables[var_name][:,:,:,:]
        else:
            var_nans = data.variables[var_name][:,:,latidxs[0]:latidxs[1],lonidxs[0]:lonidxs[1]]
    var_nans[var_nans == data.variables[var_name]._FillValue] = np.nan

    if var_dim == 3:
        var = np.nanmean(var_nans, axis=(1,2))
    elif var_dim == 4:
        var = np.nanmean(var_nans, axis=(2,3))


    scale, units = get_land_var_scale(var_name)
    if scale == -273.15:
        var = var*1.0 + scale
    else:
        var = var*scale

    data.close()

    return var

def calc_all_in_one(file_mark,file_path,var_name,var_dim,case_name,layer=None,latidxs=None,lonidxs=None):
    # day_sum  = 11020
    # if var_dim == 3:
    #     var = np.zeros([day_sum])
    # elif var_dim == 4:
    #     var = np.zeros([day_sum, layer])

    if latidxs==None or lonidxs==None:
        var = read_LIS_CABLE(file_path,var_name,var_dim) #[:,]
        print(var)
        if len(np.shape(var)) == 1:
            np.savetxt("./txt/"+file_mark+"_"+case_name+"_"+var_name+".txt",var)
        elif len(np.shape(var)) == 2:
            np.savetxt("./txt/"+file_mark+"_"+case_name+"_"+var_name+"_lvl-"+str(layer)+".txt",var[:,layer])
    else:
        var = read_LIS_CABLE(file_path,var_name,var_dim,latidxs=latidxs,lonidxs=lonidxs)
        if len(np.shape(var)) == 1:
            np.savetxt("./txt/"+file_mark+"_"+case_name+"_"+var_name+"_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+
                   "_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])+".txt",var)
        elif len(np.shape(var)) == 2:
            np.savetxt("./txt/"+file_mark+"_"+case_name+"_"+var_name+"_lvl-"+str(layer)+"_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+
                   "_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])+".txt",var[:,layer])

def plot_monthly(var_name,var_dim,year_s,year_e,layer=None):

    day_sum  = 0
    for year in np.arange(year_s+1,year_e+1):
        day_sum = day_sum + leap_year(year)
    if var_dim == 3:
        var = np.zeros([day_sum])
    elif var_dim == 4:
        var = np.zeros([day_sum, layer])
    day_s = 0

    for year in np.arange(year_s,year_e+1):
        print(year)
        if leap_year(year) == 365:
            dom = [31,28,31,30,31,30,31,31,30,31,30,31]
        else:
            dom = [31,29,31,30,31,30,31,31,30,31,30,31]
        if year == 1982:
            month = 12
            day_e = day_s + dom[month-1]
            flis  = path + 'LIS.CABLE.198212-198212.nc'
            var[day_s:day_e] = read_LIS_CABLE(flis,var_name,var_dim)
            day_s = day_e
        elif year == 2012:
            for month in np.arange(1,12):
                print(month)
                day_e = day_s + dom[month-1]
                flis = path + 'LIS.CABLE.2012'+'{:0>2}'.format(month)+'-2012'+'{:0>2}'.format(month)+'.nc'
                var[day_s:day_e] = read_LIS_CABLE(flis,var_name,var_dim)
                day_s = day_e
        else:
            for month in np.arange(1,13):
                print(month)
                day_e = day_s + dom[month-1]
                flis = path + 'LIS.CABLE.'+str(year)+'{:0>2}'.format(month)+'-'+str(year)+'{:0>2}'.format(month)+'.nc'
                var[day_s:day_e] = read_LIS_CABLE(flis,var_name,var_dim)
                day_s = day_e
        plot_time_series(var,var_name)

def plot_all_in_one(file_mark,var_name,case_names,layer=None,latidxs=None,lonidxs=None):

    # get var shape
    txt_name = "./txt/"+file_mark+"_"+case_names[0]+"_"+var_name

    if layer != None:
        txt_name = txt_name + "_lvl-"+str(layer)
    if latidxs!=None and lonidxs!=None:
        txt_name = txt_name +"_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+"_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])
    txt_name = txt_name + ".txt"

    var_tmp = np.loadtxt(txt_name,dtype='float')

    var = np.zeros((len(case_names),len(var_tmp)))
    var_tmp = None

    # read var
    for i, case_name in enumerate(case_names):
        txt_name = "./txt/"+file_mark+"_"+case_name+"_"+var_name
        if layer != None:
            txt_name = txt_name + "_lvl-"+str(layer)
        if latidxs!=None and lonidxs!=None:
            txt_name = txt_name +"_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+"_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])
        txt_name = txt_name + ".txt"

        var[i,:] = np.loadtxt(txt_name,dtype='float')

    message = case_names[0]+"_vs_"+case_names[1]
    if layer != None:
        message = message + "_lvl-"+str(layer)
    if latidxs!=None and lonidxs!=None:
        message = message + "_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+"_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])
    plot_time_series(var,var_name,case_names,message)

def plot_diff_all_in_one(file_mark,var_name,case_names,layer=None,latidxs=None,lonidxs=None):

    # read var txt file
    txt_name = var_name
    if layer != None:
        txt_name = txt_name + "_lvl-"+str(layer)
    if latidxs!=None and lonidxs!=None:
        txt_name = txt_name +"_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+"_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])

    txt_name1 = "./txt/"+file_mark+"_"+case_names[0]+"_"+txt_name+".txt"
    txt_name2 = "./txt/"+file_mark+"_"+case_names[1]+"_"+txt_name+".txt"

    var1 = np.loadtxt(txt_name1,dtype='float')
    var2 = np.loadtxt(txt_name2,dtype='float')

    # calc diff
    var_diff = var2-var1
    print(var_diff.dtype)

    message = case_names[0]+"_vs_"+case_names[1]+"_diff"
    if layer != None:
        message = message + "_lvl-"+str(layer)
    if latidxs!=None and lonidxs!=None:
        message = message + "_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+"_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])

    plot_time_series(var_diff,var_name,case_names,message)

def plot_value_diff_all_in_one(file_mark,var_name,case_names,layer=None,latidxs=None,lonidxs=None):

    # read var txt file
    txt_name = var_name
    if layer != None:
        txt_name = txt_name + "_lvl-"+str(layer)
    if latidxs!=None and lonidxs!=None:
        txt_name = txt_name +"_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+"_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])

    txt_name1 = "./txt/"+file_mark+"_"+case_names[0]+"_"+txt_name+".txt"
    txt_name2 = "./txt/"+file_mark+"_"+case_names[1]+"_"+txt_name+".txt"

    var1 = np.loadtxt(txt_name1,dtype='float')
    var2 = np.loadtxt(txt_name2,dtype='float')

    # calc diff
    var_diff = var2-var1
    print(var_diff.dtype)

    # set array FD,CTL,DIFF
    var = np.zeros((len(case_names)+1,len(var1)))
    var[0,:] = var1
    var[1,:] = var2
    var[2,:] = var_diff

    message = case_names[0]+"_vs_"+case_names[1]+"_value-diff"
    if layer != None:
        message = message + "_lvl-"+str(layer)
    if latidxs!=None and lonidxs!=None:
        message = message + "_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+"_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])
    case_names.append("Ctl-FD")

    plot_time_series(var,var_name,case_names,message)

def plot_time_series(var, var_name, case_names, message, val_min=None,val_max=None):

    fig = plt.figure(figsize=(7.2,4.5))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)

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

    colors = ['red','blue','black']#cm.Set2(np.arange(0,4))
    x = np.arange(1,len(var[0,:])+1)

    print(x)
    ax = fig.add_subplot(111)
    if len(np.shape(var)) == 1:
        ax.plot(x,var, c='red', lw=1.5, ls="-", label=case_names[1],alpha=1.)#,rasterized=True)
    elif len(np.shape(var)) == 2:
        for i in np.arange(len(var[:,0])):
            if i <  2 :
                ax.plot(x,var[i], c=colors[i], lw=1.5, ls="-", label=case_names[i], alpha=1.)#,rasterized=True)
            if i == 2 :
                ax1 = ax.twinx()
                ax1.plot(x,var[i], c=colors[i], lw=1.5, ls="-", label=case_names[i], alpha=1.)#,rasterized=True)
                # ax1.set_ylim(ranges_diff[i])
                # ax1.set_ylabel(ylabels_R[i],fontsize=14)
                # ax1.yaxis.set_tick_params(labelsize=12)

    # ax.set(xticks=x[::100], xticklabels=x[::100])
    ax.axis('tight')
    if val_min != None and val_max != None:
        ax.set_ylim(val_min,val_max)
    # ax.set_xlim(day_start,day_end)
    # ax.axvline(x=1 , ls="--")
    ax.set_ylabel(var_name)
    # ax.text(0.02, 0.95, '(b)', transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax.legend(loc='upper center', frameon=False)

    scale, units = get_land_var_scale(var_name)
    if units == None:
        plt.title(var_name, size=16)
    else:
        plt.title(var_name+' ('+units+')', size=16)

    fig.savefig("./plots/time_series/time_series_spinup_"+var_name+"_"+message , bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    plot_type = "value-diff"#"value" #"diff"
    var_3D_names =  [ "Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg","Snowf_tavg",
                       "Rainf_tavg","Evap_tavg","Qs_tavg","Qsb_tavg","VegT_tavg","AvgSurfT_tavg",
                       "Albedo_inst","SWE_inst","SnowDepth_inst","SoilWet_inst","ECanop_tavg","TVeg_tavg",
                       "FWsoil_tavg","ESoil_tavg","CanopInt_inst","SnowCover_inst","GPP_tavg","Wind_f_inst",
                       "Rainf_f_inst","Tair_f_inst", "Qair_f_inst","Psurf_f_inst","SWdown_f_inst","LWdown_f_inst"]

    var_2D_names =  ["Landmask_inst","Landcover_inst","Soiltype_inst","SandFrac_inst","ClayFrac_inst","SiltFrac_inst",
                     "SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst","Hyds_inst","Bch_inst","Sucs_inst",
                     "Elevation_inst","LAI_inst"]

    var_4D_names =  ["RelSMC_inst","SoilMoist_inst","SoilTemp_inst","SmLiqFrac_inst","SmFrozFrac_inst"]

    var_3D_basic_names = ['Evap_tavg',"ESoil_tavg","ECanop_tavg",'TVeg_tavg',"FWsoil_tavg","Qle_tavg",
                          "Qh_tavg","Qg_tavg","VegT_tavg","WaterTableD_tavg"]

    ### plot actual value ###
    if plot_type == "value":
        case_names = ['free_drain_11Jul','ctl_11Jul']
        file_name  = "LIS.CABLE.198212-201301.nc"
        file_mark  = "LIS.CABLE.198212-201301"
        latidxs    = [10,80]
        lonidxs    = [120,200]
        layer      = 6
        var_dim    = 4
        var_names  =  ["SoilMoist_inst"]

        file_paths = []
        for case_name in case_names:
            path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
            file_path  = path + file_name
            file_paths.append(file_path)

        # plot per var
        for var_name in var_names:
            for i in np.arange(len(file_paths)):
                calc_file_all_in_one(file_mark,file_paths[i],var_name,var_dim,case_names[i],layer=layer) #,latidxs=latidxs,lonidxs=lonidxs
            plot_all_in_one(file_mark,var_name,case_names)

        # plot per case
        # for case_name in case_names:
        #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
        #     file_path  = path + file_name
        #     for var_name in var_names:
        #         print(var_name)
        #         if file_type == 'monthly':
        #             year_s    = 1982
        #             year_e    = 2012
        #             plot_monthly(var_name,var_dim,year_s,year_e,layer)
        #         elif file_type == 'all_in_one':
        #             calc_all_in_one(file_path,var_name,var_dim,case_name,layer=layer,latidxs=latidxs,lonidxs=lonidxs)

    ### plot diff ###
    if plot_type == "diff":

        lat_sum  = 179
        lon_sum  = 199

        case_names = ['free_drain_11Jul','ctl_11Jul']
        file_name  = "LIS.CABLE.198212-201301.nc"
        file_mark  = "LIS.CABLE.198212-201301"
        file_paths = []
        for case_name in case_names:
            path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
            file_path  = path + file_name
            file_paths.append(file_path)

        layer     = 6
        var_dim   = 4
        var_names = ["SoilMoist_inst"]

        for var_name in var_names:
            print(var_name)
            io_opt   =  'read' #'plot' # "read"

            plot_diff_all_in_one(file_mark,file_paths,case_names,var_name,var_dim,layer)

    # plot actual value and diff
    if plot_type == "value-diff":
        case_names = ['ctl_25Jul','ctl_14Aug']
        file_name  = "LIS.CABLE.198212-201301.nc" #"LIS.CABLE.201201-201301.nc"
        file_mark  = "LIS.CABLE.198212-201301"
        latidxs    = [10,80]
        lonidxs    = [120,200]
        layer      = None
        var_dim    = 3
        var_names  = ["SoilTemp_inst", "SoilMoist_inst"] #var_3D_basic_names #["SoilTemp_inst", "SoilMoist_inst"] #

        file_paths = []
        for case_name in case_names:
            path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
            file_path  = path + file_name
            file_paths.append(file_path)

        # plot per var
        for layer in np.arange(0,6):
            for var_name in var_names:
                for i in np.arange(len(file_paths)):
                    calc_all_in_one(file_mark,file_paths[i],var_name,var_dim,case_names[i],layer=layer,latidxs=latidxs,lonidxs=lonidxs)
                plot_value_diff_all_in_one(file_mark,var_name,case_names,layer=layer,latidxs=latidxs,lonidxs=lonidxs)
