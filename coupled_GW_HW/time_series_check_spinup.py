#!/usr/bin/python

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

def plot_time_series(var,var_name,lvl):

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
    ax.plot(x, var, c=colors[0], lw=1.5, ls="-", label=var_name, alpha=1.)
    print(var)
    ax.set(xticks=x, xticklabels=x)
    ax.axis('tight')
    ax.set_ylim(0.,0.4)
    #ax.set_xlim(day_start,day_end)
    # ax.axvline(x=1 , ls="--")
    ax.set_ylabel(var_name)
    #ax.text(0.02, 0.95, '(b)', transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax.legend(loc='upper center', frameon=False)
    fig.savefig("./plots/time_series_check_spinup_"+var_name+"_"+lvl , bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    path     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hires_r7264/LIS_output/'
    var_name = "SoilMoist_inst"
    year_s   = 1982
    year_e   = 2012
    var_dim  = 4
    layer    = 6
    # lat_sum  = 179
    # lon_sum  = 199
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

    print(var)

    plot_time_series(var[:,0],var_name,"0")
    plot_time_series(var[:,1],var_name,"1")
    plot_time_series(var[:,2],var_name,"2")
    plot_time_series(var[:,3],var_name,"3")
    plot_time_series(var[:,4],var_name,"4")
    plot_time_series(var[:,5],var_name,"5")
