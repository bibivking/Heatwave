#!/usr/bin/python

'''
Plot spitial map of land diagnosis and parameters from LIS-CABLE
1. per time step
2. time period average
'''

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_range_diff
from common_utils import *

def plot_spatial_land_first_version(file_paths, wrf_path, time_s, time_e, seconds=None,loc_lat=None, loc_lon=None, message=None):

    # ==================== Read netcdf file ====================
    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    file1 = Dataset(file_paths[0], mode='r')
    Time  = nc.num2date(file1.variables['time'][:],file1.variables['time'].units,
            only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time  = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)

    file1 = Dataset(file_paths[0], mode='r')

    # "Tmax"
    T1    = file1.variables["Tair_f_inst"][:]
    tmax1 = spital_var_max(time,T1,time_s,time_e,seconds)

    # "Tmin"
    tmin1 = spital_var_min(time,T1,time_s,time_e,seconds)

    # "Qle":
    Qle1 = file1.variables["Qle_tavg"][:]
    qle1 = spital_var(time,Qle1,time_s,time_e,seconds)

    # "Qh":
    Qh1 = file1.variables["Qh_tavg"][:]
    qh1 = spital_var(time,Qh1,time_s,time_e,seconds)

    # "Rnet":
    Rnet1 = file1.variables["Qle_tavg"][:] + file1.variables["Qh_tavg"][:] + file1.variables["Qg_tavg"][:]
    rnet1 = spital_var(time,Rnet1,time_s,time_e,seconds)

    # "Lwnet_tavg":
    Lwnet1 = file1.variables["Lwnet_tavg"][:]
    lwnet1 = spital_var(time,Lwnet1,time_s,time_e,seconds)


    if len(file_paths) > 1:
        file2 = Dataset(file_paths[1], mode='r')

        # "Tmax"
        T2    = file2.variables["Tair_f_inst"][:]
        tmax2 = spital_var_max(time,T2,time_s,time_e,seconds)

        # "Tmin"
        tmin2 = spital_var_min(time,T2,time_s,time_e,seconds)

        # "Qle":
        Qle2 = file2.variables["Qle_tavg"][:]
        qle2 = spital_var(time,Qle2,time_s,time_e,seconds)

        # "Qh":
        Qh2 = file2.variables["Qh_tavg"][:]
        qh2 = spital_var(time,Qh2,time_s,time_e,seconds)

        # "Rnet":
        Rnet2 = file2.variables["Qle_tavg"][:] + file2.variables["Qh_tavg"][:] + file2.variables["Qg_tavg"][:]
        rnet2 = spital_var(time,Rnet2,time_s,time_e,seconds)

        # "Lwnet_tavg":
        Lwnet2 = file2.variables["Lwnet_tavg"][:]
        lwnet2 = spital_var(time,Lwnet2,time_s,time_e,seconds)

        tmax  = tmax2 - tmax1
        tmin  = tmin2 - tmin1
        qle   = qle2 - qle1
        qh    = qh2 - qh1
        rnet  = rnet2 - rnet1
        lwnet = lwnet2 - lwnet1
    else:
        tmax  = tmax1
        tmin  = tmin1
        qle   = qle1
        qh    = qh1
        rnet  = rnet1
        lwnet = lwnet1

    # ======================= Make plots ========================
    # Three integers (nrows, ncols, index)

    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=[9,12],sharex=True, sharey=True, squeeze=True,
                           subplot_kw={'projection': ccrs.PlateCarree()})
    plt.subplots_adjust(wspace=0.0, hspace=-0.1)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color']     = almost_black
    plt.rcParams['xtick.color']     = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']      = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor']  = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    # set the box type of sequence number
    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
    # choose colormap

    states = NaturalEarthFeature(category="cultural", scale="50m",
                                        facecolor="none",
                                        name="admin_1_states_provinces_shp")

    cmap = plt.cm.seismic #cmap = plt.cm.seismic_r

    for i in np.arange(3):
        for j in np.arange(2):
            ax[i,j].coastlines(resolution="50m",linewidth=1)
            ax[i,j].set_extent([130,155,-44,-20])
            ax[i,j].add_feature(states, linewidth=.5, edgecolor="black")

            # Add gridlines
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
            gl.xlabels_top  = False
            gl.ylabels_right= False
            gl.xlines       = False
            gl.ylines       = False
            gl.xlocator     = mticker.FixedLocator([125,130,135,140,145,150,155,160])
            gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15])
            gl.xformatter   = LONGITUDE_FORMATTER
            gl.yformatter   = LATITUDE_FORMATTER
            gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
            gl.ylabel_style = {'size':12, 'color':almost_black}
            if j == 0:
                gl.ylabels_left = True
            else:
                gl.ylabels_left = False
            if i == 2:
                gl.xlabels_bottom = True
            else:
                gl.xlabels_bottom = False

            if i == 0 and j == 0:
                clevs = np.linspace(-2,2, num=11)
                var   = tmax
                cmap  = plt.cm.seismic #cmap = plt.cm.seismic_r
                ax[i,j].text(0.02, 0.15, '(a) Tmax', transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)

            elif i == 0 and j == 1:
                clevs = np.linspace(-2,2, num=11)
                # np.linspace(-4,4, num=21)
                #[-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3.,3.5,4]
                var   = tmin
                cmap  = plt.cm.seismic
                ax[i,j].text(0.02, 0.15, '(b) Tmin', transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)

            elif i == 1 and j == 0:
                clevs = np.linspace(-80, 80, num=17)
                # np.linspace(-80, 80, num=25)
                var   = qle
                cmap  = plt.cm.seismic_r
                ax[i,j].text(0.02, 0.15, '(c) Qle', transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)

            elif i == 1 and j == 1:
                clevs = np.linspace(-80, 80, num=17)
                #np.linspace(-60, 60, num=25)
                # [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50]
                var   = qh
                cmap = plt.cm.seismic_r
                ax[i,j].text(0.02, 0.15, '(d) Qh', transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)

            elif i == 2 and j == 0:
                clevs = np.linspace(-20, 20, num=9)
                var   = rnet
                cmap  = plt.cm.seismic
                ax[i,j].text(0.02, 0.15, '(e) Rnet', transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)

            elif i == 2 and j == 1:
                clevs = np.linspace(-20, 20, num=9)
                #[-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50]
                var   = lwnet
                cmap = plt.cm.seismic
                ax[i,j].text(0.02, 0.15, '(g) LWnet', transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)

            plot = ax[i,j].contourf(lon, lat, var, levels=clevs[clevs != 0], transform=ccrs.PlateCarree(),cmap=cmap,extend='both') # [clevs != 0]
            var = None
            clevs = None

        plt.colorbar(plot, ax=ax[i], ticklocation="right", pad=0.01, orientation="vertical", aspect=20, shrink=0.88) # cax=cax,


    # plt.title(var_name, size=16)
    # cb.ax.tick_params(labelsize=10)

    plt.savefig('./plots/figures/spatial_map_'+message+'.png',dpi=300)

def plot_spatial_land(file_paths, wrf_path, time_s, time_e, seconds=None,loc_lat=None, loc_lon=None, message=None):

    # ==================== Read netcdf file ====================
    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    file1 = Dataset(file_paths[0], mode='r')
    Time  = nc.num2date(file1.variables['time'][:],file1.variables['time'].units,
            only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time  = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)

    file1 = Dataset(file_paths[0], mode='r')

    # "Tmax"
    T1    = file1.variables["Tair_f_inst"][:]
    tmax1 = spital_var_max(time,T1,time_s,time_e,seconds)

    # "Tmin"
    tmin1 = spital_var_min(time,T1,time_s,time_e,seconds)

    # "Qle":
    Qle1 = file1.variables["Qle_tavg"][:]
    qle1 = spital_var(time,Qle1,time_s,time_e,seconds)

    # "Qh":
    Qh1 = file1.variables["Qh_tavg"][:]
    qh1 = spital_var(time,Qh1,time_s,time_e,seconds)

    # "Rnet":
    Rnet1 = file1.variables["Qle_tavg"][:] + file1.variables["Qh_tavg"][:] + file1.variables["Qg_tavg"][:]
    rnet1 = spital_var(time,Rnet1,time_s,time_e,seconds)

    # "Lwnet_tavg":
    Lwnet1 = file1.variables["Lwnet_tavg"][:]
    lwnet1 = spital_var(time,Lwnet1,time_s,time_e,seconds)


    if len(file_paths) > 1:
        file2 = Dataset(file_paths[1], mode='r')

        # "Tmax"
        T2    = file2.variables["Tair_f_inst"][:]
        tmax2 = spital_var_max(time,T2,time_s,time_e,seconds)

        # "Tmin"
        tmin2 = spital_var_min(time,T2,time_s,time_e,seconds)

        # "Qle":
        Qle2 = file2.variables["Qle_tavg"][:]
        qle2 = spital_var(time,Qle2,time_s,time_e,seconds)

        # "Qh":
        Qh2 = file2.variables["Qh_tavg"][:]
        qh2 = spital_var(time,Qh2,time_s,time_e,seconds)

        # "Rnet":
        Rnet2 = file2.variables["Qle_tavg"][:] + file2.variables["Qh_tavg"][:] + file2.variables["Qg_tavg"][:]
        rnet2 = spital_var(time,Rnet2,time_s,time_e,seconds)

        # "Lwnet_tavg":
        Lwnet2 = file2.variables["Lwnet_tavg"][:]
        lwnet2 = spital_var(time,Lwnet2,time_s,time_e,seconds)

        tmax  = tmax2 - tmax1
        tmin  = tmin2 - tmin1
        qle   = qle2 - qle1
        qh    = qh2 - qh1
        rnet  = rnet2 - rnet1
        lwnet = lwnet2 - lwnet1
    else:
        tmax  = tmax1
        tmin  = tmin1
        qle   = qle1
        qh    = qh1
        rnet  = rnet1
        lwnet = lwnet1

    # ======================= Make plots ========================
    # Three integers (nrows, ncols, index)

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[9,8],sharex=True, sharey=True, squeeze=True,
                           subplot_kw={'projection': ccrs.PlateCarree()})
    plt.subplots_adjust(wspace=0.11, hspace=-0.2)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color']     = almost_black
    plt.rcParams['xtick.color']     = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']      = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor']  = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    # set the box type of sequence number
    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
    # choose colormap

    states = NaturalEarthFeature(category="cultural", scale="50m",
                                        facecolor="none",
                                        name="admin_1_states_provinces_shp")

    cmap = plt.cm.seismic #cmap = plt.cm.seismic_r

    for i in np.arange(2):
        for j in np.arange(2):
            ax[i,j].coastlines(resolution="50m",linewidth=1)
            ax[i,j].set_extent([130,155,-44,-20])
            ax[i,j].add_feature(states, linewidth=.5, edgecolor="black")

            # Add gridlines
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
            gl.xlabels_top  = False
            gl.ylabels_right= False
            gl.xlines       = False
            gl.ylines       = False
            gl.xlocator     = mticker.FixedLocator([125,130,135,140,145,150,155,160])
            gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15])
            gl.xformatter   = LONGITUDE_FORMATTER
            gl.yformatter   = LATITUDE_FORMATTER
            gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
            gl.ylabel_style = {'size':12, 'color':almost_black}
            if j == 0:
                gl.ylabels_left = True
            else:
                gl.ylabels_left = False
            if i == 1:
                gl.xlabels_bottom = True
            else:
                gl.xlabels_bottom = False

            if i == 0 and j == 0:
                clevs = np.linspace(-2,2, num=11)
                var   = tmax
                cmap  = plt.cm.seismic #cmap = plt.cm.seismic_r
                color_label= "$\mathregular{^o}$C"
                ax[i,j].text(0.02, 0.15, '(a) Tmax', transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)

            elif i == 0 and j == 1:
                clevs = np.linspace(-80, 80, num=17)
                # np.linspace(-80, 80, num=25)
                var   = qh
                cmap  = plt.cm.seismic
                color_label= "W m$\mathregular{^{-2}}$"
                ax[i,j].text(0.02, 0.15, '(b) Qh', transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)

            elif i == 1 and j == 0:
                clevs = np.linspace(-20, 20, num=9)
                var   = rnet
                cmap  = plt.cm.seismic
                color_label= "W m$\mathregular{^{-2}}$"
                ax[i,j].text(0.02, 0.15, '(c) Rnet', transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)

            elif i == 1 and j == 1:
                clevs = np.linspace(-20, 20, num=9)
                #[-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50]
                var   = lwnet
                cmap = plt.cm.seismic
                color_label= "W m$\mathregular{^{-2}}$"
                ax[i,j].text(0.02, 0.15, '(d) LWnet', transform=ax[i,j].transAxes, fontsize=14, verticalalignment='top', bbox=props)

            plot = ax[i,j].contourf(lon, lat, var, levels=clevs[clevs != 0], transform=ccrs.PlateCarree(),cmap=cmap,extend='both') # [clevs != 0]

            var = None
            clevs = None

            cbar = plt.colorbar(plot, ax=ax[i,j], ticklocation="right", pad=0.01, orientation="vertical", aspect=20, shrink=0.7) # cax=cax,
            cbar.set_label(color_label,  loc='center',size=12)
            cbar.ax.tick_params(labelsize=12)


    # plt.title(var_name, size=16)
    # cb.ax.tick_params(labelsize=10)

    plt.savefig('./plots/figures/spatial_map_'+message+'.png',dpi=300)

if __name__ == "__main__":


    # =============================== Operation ================================
    case_names = ["hw2009_3Nov", "hw2013_3Nov", "hw2019_3Nov"]

    for case_name in case_names:

        if case_name == "hw2009_3Nov":
            period     = "20090122-20090213"
            time_s = datetime(2009,1,28,0,0,0,0)
            time_e = datetime(2009,2,8,23,59,0,0)
        elif  case_name == "hw2013_3Nov":
            period     = "20121229-20130122"
            time_s = datetime(2013,1,4,0,0,0,0)
            time_e = datetime(2013,1,18,23,59,0,0)
        elif  case_name == "hw2019_3Nov":
            period     = "20190108-20190130"
            time_s = datetime(2019,1,14,0,0,0)
            time_e = datetime(2019,1,26,23,59,0,0)


        wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/ensemble_avg/wrfout_"+period+"_gw"
        cpl_land_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'

        cpl_land_file_gw  = cpl_land_file + '/LIS.CABLE.'+period+'_gw.nc'  # land output of wrf-cable run
        cpl_land_file_fd  = cpl_land_file + '/LIS.CABLE.'+period+'_fd.nc'  # land output of wrf-cable run

        file_paths        = [cpl_land_file_fd,cpl_land_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw
        seconds           = None #[18.*60.*60.,6.*60.*60.]


        if seconds == None:
            if len(file_paths) > 1:
                message = "GW-FD_"+str(time_s)+"-"+str(time_e)
            else:
                message = "GW_"+str(time_s)+"-"+str(time_e)
        else:
            if seconds[0] < seconds[1]:
                day_or_night = "Day"
            else:
                day_or_night = "Night"
            if len(file_paths) > 1:
                message = day_or_night+"_GW-FD_"+str(time_s)+"-"+str(time_e)
            else:
                message = day_or_night+"_Land_GW_"+str(time_s)+"-"+str(time_e)

        plot_spatial_land(file_paths, wrf_path, time_s, time_e, seconds=seconds, message=message)
