#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def leap_year(year):
   if (year % 4) == 0:
       return 366
   else:
       return 365

def read_LIS_CABLE(flis,var_name,var_dim):
    data = Dataset(flis, mode='r')
    var_tmp = data.variables[var_name][:]
    if var_dim == 3:
        var = np.nanmean(var_tmp, axis=(1,2))
    elif var_dim == 4:
        var = np.nanmean(var_tmp, axis=(2,3))
    data.close()
    return var

def plot_time_series(var,var_name,case_name,val_min=None,val_max=None,lvl=None):

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

    colors = cm.Set2(np.arange(0,4))
    x = np.arange(1,len(var)+1)

    print(x)
    ax = fig.add_subplot(111)
    ax.plot(x, var, c=colors[0], lw=1.5, ls="-", label=var_name, alpha=1.,rasterized=True)
    print(var)
    ax.set(xticks=x, xticklabels=x)
    ax.axis('tight')
    if val_min != None and val_max != None:
        ax.set_ylim(val_min,val_max)
    # ax.set_xlim(day_start,day_end)
    # ax.axvline(x=1 , ls="--")
    ax.set_ylabel(var_name)
    #ax.text(0.02, 0.95, '(b)', transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax.legend(loc='upper center', frameon=False)
    if lvl == None:
        fig.savefig("./plots/time_series_check_spinup_"+var_name+"_"+case_name , bbox_inches='tight', pad_inches=0.1)
    else:
        fig.savefig("./plots/time_series_check_spinup_"+var_name+"_"+lvl+"_"+case_name , bbox_inches='tight', pad_inches=0.1)

def plot_file_type_monthly(var_name,var_dim,year_s,year_e,layer=None):

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

def plot_file_type_all_in_one(var_name,var_dim,case_name,io_opt,layer=None):
    if io_opt == 'read':
        day_sum  = 11020
        if var_dim == 3:
            var = np.zeros([day_sum])
        elif var_dim == 4:
            var = np.zeros([day_sum, layer])
        flis = path + 'LIS.CABLE.198212-201301.nc'
        var[:,] = read_LIS_CABLE(flis,var_name,var_dim)
        var.dtype
        np.savetxt("./txt/"+var_name+"_"+case_name+".txt",var)
    if io_opt == 'plot':
        var = np.loadtxt("./txt/"+var_name+"_"+case_name+".txt",dtype='float')
        print(var.dtype)
        # val_min, val_max = 0. , 0.4
        plot_time_series(var,var_name,case_name)

def plot_file_type_diff_all_in_one(file_paths,case_names,var_name,var_dim,io_opt,layer):
    if io_opt == 'read':
        day_sum  = 11020
        if var_dim == 3:
            var_1 = np.zeros([day_sum])
            var_2 = np.zeros([day_sum])
        elif var_dim == 4:
            var_1 = np.zeros([day_sum, layer])
            var_2 = np.zeros([day_sum, layer])
        
        var_1[:,] = read_LIS_CABLE(file_paths[0],var_name,var_dim) 
        var_2[:,] = read_LIS_CABLE(file_paths[1],var_name,var_dim)
        var_diff  = var_2 - var_1 
        np.savetxt("./txt/"+var_name+"_"+case_names[0]+"_vs_"+case_names[1]+".txt",var_diff)
    if io_opt == 'plot':
        var = np.loadtxt("./txt/"+var_name+"_"+case_names[0]+"_vs_"+case_names[1]+".txt",dtype='float')
        case_name = case_names[0]+"_vs_"+case_names[1]
        print(var.dtype)
        # val_min, val_max = 0. , 0.4
        plot_time_series(var,var_name,case_name)

        
if __name__ == "__main__":

    ### Single Case ###

    '''
    path      = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/free_drain_hires_r7264/LIS_output/'
    case_name = "FD"
    file_type = 'all_in_one'
    layer     = 6
    var_dim   = 3
    var_names = ["Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg","Snowf_tavg",
                "Rainf_tavg","Evap_tavg","Qs_tavg","Qsb_tavg","VegT_tavg","AvgSurfT_tavg",
                "Albedo_inst","SWE_inst","SnowDepth_inst","SoilWet_inst","ECanop_tavg","TVeg_tavg",
                "FWsoil_tavg","ESoil_tavg","CanopInt_inst","SnowCover_inst","GPP_tavg","Wind_f_inst",
                "Rainf_f_inst","Tair_f_inst", "Qair_f_inst","Psurf_f_inst","SWdown_f_inst","LWdown_f_inst"]

    
    # ["Landmask_inst","Landcover_inst","Soiltype_inst","SandFrac_inst","ClayFrac_inst","SiltFrac_inst",
    # "SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst","Hyds_inst","Bch_inst","Sucs_inst",
    # "Elevation_inst","LAI_inst"]
    # ["RelSMC_inst","SoilMoist_inst","SoilTemp_inst","SmLiqFrac_inst","SmFrozFrac_inst"]
    
    # lat_sum  = 179
    # lon_sum  = 199
    for var_name in var_names:
        print(var_name)
        if file_type == 'monthly':
            year_s    = 1982
            year_e    = 2012
            plot_file_type_monthly(var_name,var_dim,year_s,year_e,layer)

        elif file_type == 'all_in_one':
            io_opt   =  'plot' # "read"
            plot_file_type_all_in_one(var_name,var_dim,case_name,io_opt,layer)
    '''


    ### Two Case ###

    case_names = ['free_drain_hires_r7264','hires_r7264']
    file_name  = "LIS.CABLE.198212-201301.nc"
    file_paths = []
    for case_name in case_names:
        path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
        file_path  = path + file_name
        file_paths.append(file_path)

    file_type = 'all_in_one'
    layer     = 6
    var_dim   = 3
    var_names = ["Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg","Snowf_tavg",
                "Rainf_tavg","Evap_tavg","Qs_tavg","Qsb_tavg","VegT_tavg","AvgSurfT_tavg",
                "Albedo_inst","SWE_inst","SnowDepth_inst","SoilWet_inst","ECanop_tavg","TVeg_tavg",
                "FWsoil_tavg","ESoil_tavg","CanopInt_inst","SnowCover_inst","GPP_tavg","Wind_f_inst",
                "Rainf_f_inst","Tair_f_inst", "Qair_f_inst","Psurf_f_inst","SWdown_f_inst","LWdown_f_inst"]

    '''
    ["Landmask_inst","Landcover_inst","Soiltype_inst","SandFrac_inst","ClayFrac_inst","SiltFrac_inst",
    "SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst","Hyds_inst","Bch_inst","Sucs_inst",
    "Elevation_inst","LAI_inst"]
    ["RelSMC_inst","SoilMoist_inst","SoilTemp_inst","SmLiqFrac_inst","SmFrozFrac_inst"]
    '''
    for var_name in var_names:
        print(var_name)
        io_opt   =  'plot' # "read"         
        plot_file_type_diff_all_in_one(file_paths,case_names,var_name,var_dim,io_opt,layer)