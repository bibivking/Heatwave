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
from cartopy.feature import NaturalEarthFeature, OCEAN
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_range_diff
from common_utils import *

def plot_spatial_fwsoil_qle_qh_tmax_qair(case_names, wrf_path, loc_lat=None, loc_lon=None, seconds=None, message=None):

    '''
    Figure 2
    '''

    # ======================= Make plots ========================
    # Three integers (nrows, ncols, index)

    fig, ax = plt.subplots(nrows=3, ncols=5, figsize=[16,10],sharex=True, sharey=True, squeeze=True,
                           subplot_kw={'projection': ccrs.PlateCarree()})
    plt.subplots_adjust(wspace=-0.44, hspace=0) # left=0.15,right=0.95,top=0.85,bottom=0.05,

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

    states= NaturalEarthFeature(category="cultural", scale="50m",
                                        facecolor="none",
                                        name="admin_1_states_provinces_shp")
    # ocean = NaturalEarthFeature('cultural', 'ocean', scale='50m',
    #                         edgecolor='none', facecolor="lightgray")

    # ======================= Set colormap =======================
    cmap       = plt.cm.seismic
    blue2white = truncate_colormap(cmap, minval=0., maxval=0.5)
    white2red  = truncate_colormap(cmap, minval=0.5, maxval=1.)
    cmap       = plt.cm.seismic_r
    red2white  = truncate_colormap(cmap, minval=0., maxval=0.5)
    white2blue = truncate_colormap(cmap, minval=0.5, maxval=1.)

    # ======================= Read WRF file =======================
    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    texts = [ "(a)","(b)","(c)","(d)","(e)",
              "(f)","(g)","(h)","(i)","(j)",
              "(k)","(l)","(m)","(n)","(o)"]

    label_x = ["Δ$β$","ΔQ$\mathregular{_e}$","ΔQ$\mathregular{_h}$","ΔT$\mathregular{_{max}}$","Δq"]
    label_y = ["2009","2013","2019"]
    loc_y   = [0.55,0.5,0.45]

    cnt     = 0

    # ==================== Set up files ====================
    for i, case_name in enumerate(case_names):

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

        cpl_land_file     = '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'

        cpl_land_file_gw  = cpl_land_file + '/LIS.CABLE.'+period+'_gw.nc'  # land output of wrf-cable run
        cpl_land_file_fd  = cpl_land_file + '/LIS.CABLE.'+period+'_fd.nc'  # land output of wrf-cable run

        file_paths        = [cpl_land_file_fd,cpl_land_file_gw]

        # ==================== Read netcdf file ====================
        # Open the NetCDF4 file (add a directory path if necessary) for reading:
        file1 = Dataset(file_paths[0], mode='r')
        Time  = nc.num2date(file1.variables['time'][:],file1.variables['time'].units,
                only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        time  = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)

        # "Tmax1"
        T1    = file1.variables["Tair_f_inst"][:]
        tmax1 = spital_var_max(time,T1,time_s,time_e,seconds)

        # "Qle1"
        Qle1  = file1.variables["Qle_tavg"][:]
        qle1  = spital_var(time,Qle1,time_s,time_e,seconds)

        # "FW1"
        FW1    = file1.variables["FWsoil_tavg"][:]
        fw1    = spital_var(time,FW1,time_s,time_e,seconds)

        # "Qh1"
        Qh1   = file1.variables["Qh_tavg"][:]
        qh1   = spital_var(time,Qh1,time_s,time_e,seconds)

        # "Q1"
        Q1    = file1.variables["Qair_f_inst"][:]
        q1    = spital_var(time,Q1,time_s,time_e,seconds)

        if len(file_paths) > 1:
            file2 = Dataset(file_paths[1], mode='r')

            # "Tmax1"
            T2    = file2.variables["Tair_f_inst"][:]
            tmax2 = spital_var_max(time,T2,time_s,time_e,seconds)

            # "Qle2"
            Qle2  = file2.variables["Qle_tavg"][:]
            qle2  = spital_var(time,Qle2,time_s,time_e,seconds)

            # "FW2"
            FW2   = file2.variables["FWsoil_tavg"][:]
            fw2   = spital_var(time,FW2,time_s,time_e,seconds)

            # "Qh2"
            Qh2   = file2.variables["Qh_tavg"][:]
            qh2   = spital_var(time,Qh2,time_s,time_e,seconds)

            # "Q2"
            Q2    = file2.variables["Qair_f_inst"][:]
            q2    = spital_var(time,Q2,time_s,time_e,seconds)

            tmax  = tmax2 - tmax1
            qle   = qle2 - qle1
            fw    = fw2 - fw1
            qh    = qh2 - qh1
            q     = q2 - q1

        else:
            qh    = qh1
            q     = q1
            tmax  = tmax1
            qle   = qle1
            fw    = fw1

        # ==================== Start to plot ====================
        for j in np.arange(5):

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
                gl.ylabels_left   = True
            else:
                gl.ylabels_left   = False
            if i == 2:
                gl.xlabels_bottom = True
            else:
                gl.xlabels_bottom = False

            # set y label
            if j == 0:
                ax[i,j].text(-0.23, loc_y[i], label_y[i], va='bottom', ha='center',
                              rotation='vertical', rotation_mode='anchor',
                              transform=ax[i,j].transAxes)

        # left - FWsoil
        clevs1  = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.]
        plot1   = ax[i,0].contourf(lon, lat, fw, levels=clevs1, transform=ccrs.PlateCarree(),cmap=white2blue,extend='both')
        ax[i,0].text(0.02, 0.15, texts[cnt], transform=ax[i,0].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,0].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # middle left - Qle
        clevs2   = [ 10, 20, 40, 60, 80, 100, 120, 140]
        plot2    = ax[i,1].contourf(lon, lat, qle, levels=clevs2, transform=ccrs.PlateCarree(),cmap=white2blue,extend='both')
        ax[i,1].text(0.02, 0.15, texts[cnt+1], transform=ax[i,1].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,1].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # middle - Qh
        clevs3  = [ -140, -120, -100, -80, -60, -40, -20, -10]
        plot3   = ax[i,2].contourf(lon, lat, qh, levels=clevs3, transform=ccrs.PlateCarree(),cmap=blue2white,extend='both')
        ax[i,2].text(0.02, 0.15, texts[cnt+2], transform=ax[i,2].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,2].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # middle right - Tmax
        tmax     = np.where(np.isnan(fw), np.nan, tmax)
        clevs4   = [-3.,-2.5,-2.,-1.5,-1.,-0.5, -0.1]
        plot4    = ax[i,3].contourf(lon, lat, tmax, levels=clevs4, transform=ccrs.PlateCarree(),cmap=blue2white,extend='both') #
        ax[i,3].text(0.02, 0.15, texts[cnt+3], transform=ax[i,3].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,3].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # right - Qair
        clevs5  = [0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.]
        plot5   = ax[i,4].contourf(lon, lat, q*1000., levels=clevs5, transform=ccrs.PlateCarree(),cmap=white2blue,extend='both')
        ax[i,4].text(0.02, 0.15, texts[cnt+4], transform=ax[i,4].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,4].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")
        clevs   = None

        # set top x label
        if i == 0:
            # ax[i,0].xaxis.set_label_position('top')
            # ax[i,1].xaxis.set_label_position('top')
            # ax[i,2].xaxis.set_label_position('top')
            # ax[i,0].set_xlabel(label_x[0],labelpad=-0.1)#, fontsize=12)
            # ax[i,1].set_xlabel(label_x[1],labelpad=-0.1)#, fontsize=12)
            # ax[i,2].set_xlabel(label_x[2],labelpad=-0.1)#, fontsize=12)

            ax[i,0].set_title(label_x[0])#, fontsize=12)
            ax[i,1].set_title(label_x[1])#, fontsize=12)
            ax[i,2].set_title(label_x[2])#, fontsize=12)
            ax[i,3].set_title(label_x[3])#, fontsize=12)
            ax[i,4].set_title(label_x[4])#, fontsize=12)


        # set bottom colorbar
        if i == 2:

            # left - Fwsoil
            cbar = plt.colorbar(plot1, ax=ax[:,0], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.49)
            color_label= "-"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=10, rotation=45)

            # middle left - Qle
            cbar = plt.colorbar(plot2, ax=ax[:,1], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.49)
            color_label= "W m$\mathregular{^{-2}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=10, rotation=45)

            # middle - Qh
            cbar = plt.colorbar(plot3, ax=ax[:,2], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.49)
            color_label= "W m$\mathregular{^{-2}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=10, rotation=45)

            # middle right - Tmax
            cbar = plt.colorbar(plot4, ax=ax[:,3], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.49)
            color_label= "$\mathregular{^o}$C"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=10, rotation=45)

            # right - q
            cbar = plt.colorbar(plot5, ax=ax[:,4], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.49)
            color_label= "g kg$\mathregular{^{-1}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=10, rotation=45)

        cnt = cnt + 5
    # plt.tight_layout(pad=0.01)
    plt.savefig('./plots/Fig2_spatial_map_fwsoil_Qle_Qh_Tmax_q_2009_2013_2019.png',dpi=300)
    #, bbox_inches='tight', pad_inches=0.1)

def plot_spatial_Rnet_LW_SW(case_names, loc_lat=None, loc_lon=None, seconds=None, message=None):

    '''
    old Figure S3, now Figure 3
    '''

    # ======================= Make plots ========================
    # Three integers (nrows, ncols, index)

    fig, ax = plt.subplots(nrows=3, ncols=5, figsize=[16,10],sharex=True, sharey=True, squeeze=True,
                           subplot_kw={'projection': ccrs.PlateCarree()})
    plt.subplots_adjust(wspace=-0.44, hspace=0) # left=0.15,right=0.95,top=0.85,bottom=0.05,

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

    states= NaturalEarthFeature(category="cultural", scale="50m",
                                        facecolor="none",
                                        name="admin_1_states_provinces_shp")
    # ocean = NaturalEarthFeature('cultural', 'ocean', scale='50m',
    #                             edgecolor='none', facecolor="lightgray")

    # ======================= Set colormap =======================
    cmap       = plt.cm.seismic

    # ======================= Read WRF file =======================
    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    texts   = [ "(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)","(o)"]
    label_x = ["ΔR$\mathregular{_{net}}$","ΔLW$\mathregular{_{net}}$","ΔLW$\mathregular{_{up}}$",
               "ΔLW$\mathregular{_{dn}}$","ΔSW$\mathregular{_{net}}$"]
    label_y = ["2009","2013","2019"]
    loc_y   = [0.53,0.5,0.48]
    cnt     = 0

    # ==================== Set up files ====================
    for i, case_name in enumerate(case_names):

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

        cpl_path          = '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'

        cpl_land_file_gw  = cpl_path + '/LIS.CABLE.'+period+'_gw.nc'  # land output of wrf-cable run
        cpl_land_file_fd  = cpl_path + '/LIS.CABLE.'+period+'_fd.nc'  # land output of wrf-cable run

        cpl_atmo_file_gw  = cpl_path + '/wrfout_'+period+'_gw'  # land output of wrf-cable run
        cpl_atmo_file_fd  = cpl_path + '/wrfout_'+period+'_fd'  # land output of wrf-cable run

        land_paths        = [cpl_land_file_fd,cpl_land_file_gw]
        atmo_paths        = [cpl_atmo_file_fd,cpl_atmo_file_gw]

        # ==================== Read netcdf file ====================
        # Open the NetCDF4 file (add a directory path if necessary) for reading:
        land1 = Dataset(land_paths[0], mode='r')
        Time  = nc.num2date(land1.variables['time'][:],land1.variables['time'].units,
                only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        time_land = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)

        # "Rnet":
        Rnet1 = land1.variables["Qle_tavg"][:] + land1.variables["Qh_tavg"][:] + land1.variables["Qg_tavg"][:]
        rnet1 = spital_var(time_land,Rnet1,time_s,time_e,seconds)

        # "Lwnet_tavg":
        Lwnet1 = land1.variables["Lwnet_tavg"][:]
        lwnet1 = spital_var(time_land,Lwnet1,time_s,time_e,seconds)

        # "Swnet_tavg":
        Swnet1 = land1.variables["Swnet_tavg"][:]
        swnet1 = spital_var(time_land,Swnet1,time_s,time_e,seconds)

        atmo1    = Dataset(atmo_paths[0])
        ntime    = len(atmo1.variables['Times'][:,0])
        time_tmp = []

        # Atmosphere
        encoding    = 'utf-8' # Times in WRF output is btype, convert to string
        for t in np.arange(ntime):
            time_temp = datetime.strptime(str(atmo1.variables['Times'][t,:], encoding),'%Y-%m-%d_%H:%M:%S')
            time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1))

        time_atmo = np.array(time_tmp)

        # "LWup"
        LWup1  = read_wrf_surf_var(atmo_paths[0],"LWUPB")
        lwup1  = spital_var(time_atmo,LWup1,time_s,time_e,seconds)

        # "LWdn"
        LWdn1  = read_wrf_surf_var(atmo_paths[0],"LWDNB")
        lwdn1  = spital_var(time_atmo,LWdn1,time_s,time_e,seconds)

        if len(land_paths) > 1:

            land2 = Dataset(land_paths[1], mode='r')
            Time  = nc.num2date(land2.variables['time'][:],land2.variables['time'].units,
                    only_use_cftime_datetimes=False, only_use_python_datetimes=True)
            time_land = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)

            # "Rnet":
            Rnet2 = land2.variables["Qle_tavg"][:] + land2.variables["Qh_tavg"][:] + land2.variables["Qg_tavg"][:]
            rnet2 = spital_var(time_land,Rnet2,time_s,time_e,seconds)

            # "Lwnet_tavg":
            Lwnet2 = land2.variables["Lwnet_tavg"][:]
            lwnet2 = spital_var(time_land,Lwnet2,time_s,time_e,seconds)

            # "Swnet_tavg":
            Swnet2 = land2.variables["Swnet_tavg"][:]
            swnet2 = spital_var(time_land,Swnet2,time_s,time_e,seconds)

            atmo2    = Dataset(atmo_paths[1])
            ntime    = len(atmo2.variables['Times'][:,0])
            time_tmp = []

            for t in np.arange(ntime):
                time_temp = datetime.strptime(str(atmo2.variables['Times'][t,:], encoding),'%Y-%m-%d_%H:%M:%S')
                time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1))

            time_atmo = np.array(time_tmp)

            # "LWup"
            LWup2  = read_wrf_surf_var(atmo_paths[1],"LWUPB")
            lwup2  = spital_var(time_atmo,LWup2,time_s,time_e,seconds)

            # "LWdn"
            LWdn2  = read_wrf_surf_var(atmo_paths[1],"LWDNB")
            lwdn2  = spital_var(time_atmo,LWdn2,time_s,time_e,seconds)

            rnet  = rnet2 - rnet1
            lwnet = lwnet2 - lwnet1
            swnet = swnet2 - swnet1
            lwup  = lwup2 - lwup1
            lwdn  = lwdn2 - lwdn1
        else:
            rnet  = rnet1
            lwnet = lwnet1
            swnet = swnet1
            lwup  = lwup1
            lwdn  = lwdn1

        # ==================== Start to plot ====================
        for j in np.arange(5):

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
                gl.ylabels_left   = True
            else:
                gl.ylabels_left   = False
            if i == 2:
                gl.xlabels_bottom = True
            else:
                gl.xlabels_bottom = False

            # set y label
            if j == 0:
                ax[i,j].text(-0.23, loc_y[i], label_y[i], va='bottom', ha='center',
                              rotation='vertical', rotation_mode='anchor',
                              transform=ax[i,j].transAxes)


        clevs          = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
        clevs_label    = ["-22","-18","-14","-10","-6","-2","2","6","10","14","18","22"]
        # clevs    = [-25,-20,-15,-10,-5,5,10,15,20,25]

        # left - Rnet
        plot1    = ax[i,0].contourf(lon, lat, rnet, levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
        ax[i,0].text(0.02, 0.15, texts[cnt], transform=ax[i,0].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,0].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # left middle - LWnet
        plot2    = ax[i,1].contourf(lon, lat, lwnet, levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
        ax[i,1].text(0.02, 0.15, texts[cnt+1], transform=ax[i,1].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,1].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # middle - LWup
        lwup = np.where(np.isnan(rnet), np.nan, lwup)
        plot3   = ax[i,2].contourf(lon, lat, lwup, levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
        ax[i,2].text(0.02, 0.15, texts[cnt+2], transform=ax[i,2].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,2].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # right middle - LWdn
        lwdn = np.where(np.isnan(rnet), np.nan, lwdn)
        plot4   = ax[i,3].contourf(lon, lat, lwdn, levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
        ax[i,3].text(0.02, 0.15, texts[cnt+3], transform=ax[i,3].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,3].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # right - SWnet
        plot5   = ax[i,4].contourf(lon, lat, swnet, levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
        ax[i,4].text(0.02, 0.15, texts[cnt+4], transform=ax[i,4].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,4].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        if i == 2 and j == 4:

            cbar = plt.colorbar(plot1, ax=ax, ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=40, shrink=0.8) # cax=cax,
            cbar.set_ticks(clevs)
            cbar.set_ticklabels(clevs_label) # cax=cax,
            color_label= "W m$\mathregular{^{-2}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=12)

        # set top x label
        if i == 0:
            ax[i,0].set_title(label_x[0])#,labelpad=-0.1)#, fontsize=12)
            ax[i,1].set_title(label_x[1])#,labelpad=-0.1)#, fontsize=12)
            ax[i,2].set_title(label_x[2])#,labelpad=-0.1)#, fontsize=12)
            ax[i,3].set_title(label_x[3])#,labelpad=-0.1)#, fontsize=12)
            ax[i,4].set_title(label_x[4])#,labelpad=-0.1)#, fontsize=12)

        cnt = cnt + 5
        # plt.title(var_name, size=16)
        # cb.ax[i,j].tick_params(labelsize=10)
    # plt.tight_layout(pad=0.01)
    plt.savefig('./plots/Fig3_spatial_map_plot_spatial_Rnet_LW_SW_2009_2013_2019.png',dpi=300)
    #, bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    # ======================================
    # Decks for plot_spatial_tmax_qle_fwsoil
    # ======================================

    case_names = ["hw2009_3Nov", "hw2013_3Nov", "hw2019_3Nov"]
    wrf_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg/wrfout_20090122-20090213_gw"
    # plot_spatial_fwsoil_qle_tmax(case_names, wrf_path)
    plot_spatial_Rnet_LW_SW(case_names, wrf_path)
    # plot_spatial_qh_qair(case_names, wrf_path)
    # plot_spatial_fwsoil_qle_qh_tmax_qair(case_names, wrf_path)
