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

def get_var_scale(var_name):

    '''
    Convert units
    '''

    var_s2d        = ["Rainf_tavg","Evap_tavg","ECanop_tavg","TVeg_tavg","ESoil_tavg","Qs_tavg","Qsb_tavg"]
    var_umol_s2g_d = "GPP_tavg"

    s2d            = 3600*24. # s-1 to d-1
    mmd2wm2        = 28.94  # 1 mm/day = 28.94 W/m2 
    umol_s2g_d     = 0.000001*12*s2d # umol s-1 to g d-1


    if var_name in var_s2d:
        scale = s2d*mmd2wm2
        units = "W/m2" #"mm d-1"
    elif var_name in var_umol_s2g_d:
        scale = umol_s2g_d
        units = "gC m-2 d-1"
    else:
        scale = 1.
        units = None

    return (scale, units)       

def read_LIS_CABLE(flis,var_name,var_dim,latidxs=None,lonidxs=None):

    data = Dataset(flis, mode='r')

    # Replace _FillValues with NaNs:
    if var_dim == 3:
        var_nans = data.variables[var_name][:,latidxs[0]:latidxs[1],lonidxs[0]:lonidxs[1]]
    elif var_dim == 4:
        var_nans = data.variables[var_name][:,:,latidxs[0]:latidxs[1],lonidxs[0]:lonidxs[1]]
    var_nans[var_nans == data.variables[var_name]._FillValue] = np.nan

    if var_dim == 3:
        var = np.nanmean(var_nans, axis=(1,2))
    elif var_dim == 4:
        var = np.nanmean(var_nans, axis=(2,3))
    data.close()
    print(var[5000])
    return var

def plot_time_series(var,var_name,case_name, units= None, val_min=None,val_max=None,lvl=None,latidxs=None,lonidxs=None):

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
    x = np.arange(1,len(var[0,:])+1)

    print(x)
    ax = fig.add_subplot(111)
    ax.plot(x,var[0], c='red', lw=1.5, ls="-", label=case_name[0], alpha=1.)#,rasterized=True)
    ax.plot(x,var[1], c='blue', lw=1.5, ls="-", label=case_name[1], alpha=1.)#,rasterized=True)
    print(var)
    # ax.set(xticks=x[::100], xticklabels=x[::100])
    ax.axis('tight')
    if val_min != None and val_max != None:
        ax.set_ylim(val_min,val_max)
    # ax.set_xlim(day_start,day_end)
    # ax.axvline(x=1 , ls="--")
    ax.set_ylabel(var_name)
    # ax.text(0.02, 0.95, '(b)', transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax.legend(loc='upper center', frameon=False)
    if len(case_name) == 1:
        case_name_string = case_name
    else:
        case_name_string = case_name[0]+"_vs_"+case_name[1]
    if units == None:
        plt.title(var_name, size=16)
    else:
        plt.title(var_name+' ('+units+')', size=16)
    if lvl == None:
        if latidxs==None or lonidxs==None:
            fig.savefig("./plots/time_series_check_spinup_"+var_name+"_"+case_name_string , bbox_inches='tight', pad_inches=0.1)
        else:
            fig.savefig("./plots/time_series_check_spinup_"+var_name+"_"+case_name_string+"_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+
                    "_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1]) , bbox_inches='tight', pad_inches=0.1)
    else:
        if latidxs==None or lonidxs==None:
            fig.savefig("./plots/time_series_check_spinup_"+var_name+"_"+lvl+"_"+case_name_string , bbox_inches='tight', pad_inches=0.1)
        else:
            fig.savefig("./plots/time_series_check_spinup_"+var_name+"_"+lvl+"_"+case_name_string+"_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+
                    "_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1]) , bbox_inches='tight', pad_inches=0.1)

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

def calc_file_type_all_in_one(file_path,var_name,var_dim,case_name,layer=None,latidxs=None,lonidxs=None):
    day_sum  = 11020
    if var_dim == 3:
        var = np.zeros([day_sum])
    elif var_dim == 4:
        var = np.zeros([day_sum, layer])
    if latidxs==None or lonidxs==None:
        var[:,] = read_LIS_CABLE(file_path,var_name,var_dim)
        np.savetxt("./txt/"+var_name+"_"+case_name+".txt",var)
    else:
        var[:,] = read_LIS_CABLE(file_path,var_name,var_dim,latidxs=latidxs,lonidxs=lonidxs)
        np.savetxt("./txt/"+var_name+"_"+case_name+"_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+
                   "_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])+".txt",var)
        
def plot_file_type_all_in_one(var_name,case_names,latidxs=None,lonidxs=None):
    if latidxs==None or lonidxs==None:
        var_tmp = np.loadtxt("./txt/"+var_name+"_"+case_names[0]+".txt",dtype='float')
    else:
        var_tmp = np.loadtxt("./txt/"+var_name+"_"+case_names[0]+"_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+
                  "_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])+".txt",dtype='float')

    var = np.zeros((len(case_names),len(var_tmp)))
    scale, units = get_var_scale(var_name)
    for i, case_name in enumerate(case_names):
        if latidxs==None or lonidxs==None:
            var[i,:] = np.loadtxt("./txt/"+var_name+"_"+case_name+".txt",dtype='float')*scale
        else:
            var[i,:] = np.loadtxt("./txt/"+var_name+"_"+case_name+"_lat-"+str(latidxs[0])+"-"+str(latidxs[1])+
                    "_lon-"+str(lonidxs[0])+"-"+str(lonidxs[1])+".txt",dtype='float')*scale
    print(var.shape)
        # val_min, val_max = 0. , 0.4
    plot_time_series(var,var_name,case_names,units,latidxs=latidxs,lonidxs=lonidxs)

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
    case_names = ['free_drain_hires_r7264','hires_r7264']
    file_name  = "LIS.CABLE.198212-201301.nc"
    file_type  = 'all_in_one'
    latidxs    = [129,130]
    lonidxs    = [109,110]
    layer      = 6
    var_dim    = 3
    var_names  = ['Evap_tavg',"ESoil_tavg","ECanop_tavg",'TVeg_tavg',"FWsoil_tavg","Qle_tavg","Qh_tavg","Qg_tavg"]
    '''
    ["Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg","Snowf_tavg",
                 "Rainf_tavg","Evap_tavg","Qs_tavg","Qsb_tavg","VegT_tavg","AvgSurfT_tavg",
                 "Albedo_inst","SWE_inst","SnowDepth_inst","SoilWet_inst","ECanop_tavg","TVeg_tavg",
                 "FWsoil_tavg","ESoil_tavg","CanopInt_inst","SnowCover_inst","GPP_tavg","Wind_f_inst",
                 "Rainf_f_inst","Tair_f_inst", "Qair_f_inst","Psurf_f_inst","SWdown_f_inst","LWdown_f_inst"]
    '''
    # for case_name in case_names:
    #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    #     file_path  = path + file_name
    #     for var_name in var_names:
    #         print(var_name)
    #         if file_type == 'monthly':
    #             year_s    = 1982
    #             year_e    = 2012
    #             plot_file_type_monthly(var_name,var_dim,year_s,year_e,layer)
    #         elif file_type == 'all_in_one':
    #             calc_file_type_all_in_one(file_path,var_name,var_dim,case_name,layer=layer,latidxs=latidxs,lonidxs=lonidxs)
    for var_name in var_names:
        if file_type == 'all_in_one':
            plot_file_type_all_in_one(var_name,case_names,latidxs=latidxs,lonidxs=lonidxs)
    
    # ["Landmask_inst","Landcover_inst","Soiltype_inst","SandFrac_inst","ClayFrac_inst","SiltFrac_inst",
    # "SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst","Hyds_inst","Bch_inst","Sucs_inst",
    # "Elevation_inst","LAI_inst"]
    # ["RelSMC_inst","SoilMoist_inst","SoilTemp_inst","SmLiqFrac_inst","SmFrozFrac_inst"]

    # lat_sum  = 179
    # lon_sum  = 199

    # ### Two Case ###

    # case_names = ['free_drain_hires_r7264','hires_r7264']
    # file_name  = "LIS.CABLE.198212-201301.nc"
    # file_paths = []
    # for case_name in case_names:
    #     path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    #     file_path  = path + file_name
    #     file_paths.append(file_path)

    # file_type = 'all_in_one'
    # layer     = 6
    # var_dim   = 3
    # var_names = ["Swnet_tavg"]

    # '''
    # "Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg","Snowf_tavg",
    # "Rainf_tavg","Evap_tavg","Qs_tavg","Qsb_tavg","VegT_tavg","AvgSurfT_tavg",
    # "Albedo_inst","SWE_inst","SnowDepth_inst","SoilWet_inst","ECanop_tavg","TVeg_tavg",
    # "FWsoil_tavg","ESoil_tavg","CanopInt_inst","SnowCover_inst","GPP_tavg","Wind_f_inst",
    # "Rainf_f_inst","Tair_f_inst", "Qair_f_inst","Psurf_f_inst","SWdown_f_inst","LWdown_f_inst"

    # ["Landmask_inst","Landcover_inst","Soiltype_inst","SandFrac_inst","ClayFrac_inst","SiltFrac_inst",
    # "SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst","Hyds_inst","Bch_inst","Sucs_inst",
    # "Elevation_inst","LAI_inst"]
    # ["RelSMC_inst","SoilMoist_inst","SoilTemp_inst","SmLiqFrac_inst","SmFrozFrac_inst"]
    # '''
    # for var_name in var_names:
    #     print(var_name)
    #     io_opt   =  'read' #'plot' # "read"
    #     plot_file_type_diff_all_in_one(file_paths,case_names,var_name,var_dim,io_opt,layer)
