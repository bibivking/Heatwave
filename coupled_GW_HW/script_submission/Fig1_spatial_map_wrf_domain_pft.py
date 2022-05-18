#!/usr/bin/env python
"""
Produce the Australia map with the label of EucFACE site - Figure 1
"""

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"


#!/usr/bin/python

import sys
import cartopy
import numpy as np
from netCDF4 import Dataset,num2date
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature
import scipy.ndimage as ndimage
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                        cartopy_ylim, latlon_coords, ALL_TIMES)
from common_utils import *

def spatial_map_wrf_domain(file_paths, wrf_path):

    # fig, ax = plt.subplots(nrows=3, ncols=1, ,sharex=True, sharey=True, squeeze=True)
    fig = plt.figure(figsize=[15,8],constrained_layout=True)
    plt.subplots_adjust(wspace=-0.65, hspace=0.07)
    gs  = fig.add_gridspec(nrows=2,ncols=3)

    ax1 = fig.add_subplot(gs[0,0:2],projection=ccrs.PlateCarree())
    ax2 = fig.add_subplot(gs[0,2],projection=ccrs.PlateCarree())
    ax3 = fig.add_subplot(gs[1,0],projection=ccrs.PlateCarree(),sharey=ax2)
    ax4 = fig.add_subplot(gs[1,1],projection=ccrs.PlateCarree(),sharey=ax2)
    ax5 = fig.add_subplot(gs[1,2],projection=ccrs.PlateCarree(),sharey=ax2,sharex=ax2)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black                    = '#262626'
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

    # =============== make domian ===============
    # Open the NetCDF file
    encoding = 'utf-8' # Times in WRF output is btype, convert to string

    file  = Dataset(file_paths[0])
    wrf   = Dataset(wrf_path)

    iveg  = file.variables['Landcover_inst'][0,:,:]
    lndmsk= file.variables['Landmask_inst'][0,:,:]
    print(np.unique(iveg))

    # to get lat and lon
    p     = getvar(wrf, "pressure")

    # Get the lat/lon coordinates
    lats, lons = latlon_coords(p)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(p)

    # ------- Plot Advection Colorfill -------
    # iveg in SE AUS : [2.0 5.0 6.0 9.0 14.0 ]
    #                  ["BEF","shrub","grass","crop","barren"]
    iveg = np.where(iveg == 2, 1, iveg)
    iveg = np.where(iveg == 5, 2, iveg)
    iveg = np.where(iveg == 6, 3, iveg)
    iveg = np.where(iveg == 9, 4, iveg)
    iveg = np.where(iveg == 14,5, iveg)
    iveg = np.where(iveg == -9999,0, iveg)
    np.savetxt("check",iveg)
    # =============== WTD ===============
    # Open the NetCDF file

    # 2009
    time_s    = datetime(2009,1,28,0,0,0,0)
    time_e    = datetime(2009,2,8,23,59,0,0)
    time, WTD = read_var(file_paths[0], 'WaterTableD_tavg', lat_name="lat", lon_name="lon")
    wtd_2009  = spital_var(time, WTD, time_s, time_e)
    time      = None
    WTD       = None

    # 2013
    time_s    = datetime(2013,1,4,0,0,0,0)
    time_e    = datetime(2013,1,18,23,59,0,0)
    time, WTD = read_var(file_paths[1], 'WaterTableD_tavg', lat_name="lat", lon_name="lon")
    wtd_2013  = spital_var(time, WTD, time_s, time_e)
    time      = None
    WTD       = None

    # 2019
    time_s    = datetime(2019,1,14,0,0,0)
    time_e    = datetime(2019,1,26,23,59,0,0)
    time, WTD = read_var(file_paths[2], 'WaterTableD_tavg', lat_name="lat", lon_name="lon")
    wtd_2019  = spital_var(time, WTD, time_s, time_e)
    time      = None
    WTD       = None

    lai_file  = Dataset(file_paths[3], mode='r')
    lai_jan   = lai_file.variables['LAI'][0,:,:]
    print("lai_jan",lai_jan)

    # mask out sea pixels
    wtd_2009       = np.where(lndmsk == 1., wtd_2009, np.nan)
    wtd_2013       = np.where(lndmsk == 1., wtd_2013, np.nan)
    wtd_2019       = np.where(lndmsk == 1., wtd_2019, np.nan)
    print("wtd_2019",wtd_2019)
    lai_jan        = np.where(lndmsk == 1., lai_jan,  np.nan)
    print("lai_jan",lai_jan)

    # ================ Plotting PFT ================
    # Set the GeoAxes to the projection used by WRF
    ax1.set_extent([109,161,-43,-9])
    ax1.coastlines(resolution="50m",linewidth=1)
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                         facecolor='lightgray',
                                         name="admin_1_states_provinces_shp")
    ax1.add_feature(states, linewidth=.5, zorder=0, edgecolor='lightgray', facecolor='lightgray') # '#fff9e6'
    ax1.add_feature(cartopy.feature.LAND, edgecolor='lightgray', facecolor='lightgray')
    ax1.add_feature(cartopy.feature.OCEAN, edgecolor='lightgray', facecolor='lightgray')
    ax1.add_feature(cartopy.feature.COASTLINE, lw=0.5)

    # Add gridlines
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color=almost_black, linestyle='--')
    gl.xlabels_top  = False
    gl.ylabels_right= False
    gl.xlines       = False
    gl.ylines       = False
    gl.xlocator     = mticker.FixedLocator([110,120,130,140,150,160])
    gl.ylocator     = mticker.FixedLocator([-40,-30,-20,-10])
    gl.xformatter   = LONGITUDE_FORMATTER
    gl.yformatter   = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':almost_black}#,'rotation': 90}
    gl.ylabel_style = {'size':10, 'color':almost_black}

    cint = [-1,0.1,1.1,2.1,3.1,4.1,5.1]
    almost_black_rgb = np.array([38./255, 38./255, 38./255, 1])

    cmap = plt.cm.Paired(np.arange(12))
            # ["BEF","shrub","grass","crop","barren"]
    cmap = [cmap[0],cmap[3],cmap[2],cmap[6],cmap[7],cmap[11]] #

    cf   = ax1.contourf(to_np(lons), to_np(lats), iveg, levels=cint, colors=cmap,
                       transform=ccrs.PlateCarree(),rasterized=True)
    cbar = plt.colorbar(cf, ax=ax1, orientation="horizontal", pad=.1,  aspect=35, extend='neither', shrink=0.4)
    cbar.set_ticks([])

    # add line
    x_values = [139,152]  # gather x-values
    y_values = [-36, -36] # gather y-values
    line = ax1.plot(x_values, y_values,c=almost_black, alpha = 1, linestyle="--", transform=ccrs.PlateCarree())

    # add square
    ax1.add_patch(Polygon([[148., -36.5], [149., -36.5], [149., -35.5], [148., -35.5]],
                 closed=True,color=almost_black, fill=False))
    ax1.text(-0.05, -0.2, "Water", transform=ax1.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    ax1.text(0.15, -0.2, "BEF", transform=ax1.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    ax1.text(0.33, -0.2, "Shrub", transform=ax1.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    ax1.text(0.53, -0.2, "Grass", transform=ax1.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    ax1.text(0.735, -0.2, "Crop", transform=ax1.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    ax1.text(0.92, -0.2, "Barren\n  land", transform=ax1.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    ax1.text(0.02, 0.1, "(a)", transform=ax1.transAxes, fontsize=12, verticalalignment='top', bbox=props)


    # ================ Plotting LAI ================
    # ax2 = plt.axes(projection=ccrs.PlateCarree())
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                        facecolor="none",
                                        name="admin_1_states_provinces_shp")
    ax2.set_extent([130,155,-44,-20])
    ax2.coastlines(resolution="50m",linewidth=1)
    ax2.add_feature(states, linewidth=.5, edgecolor="black")
    ax2.add_feature(cartopy.feature.OCEAN, edgecolor='lightgray', facecolor='lightgray')

    # Add gridlines
    gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
    gl.xlabels_top  = False
    gl.ylabels_right= False
    gl.xlines       = False
    gl.ylines       = False
    gl.xlocator     = mticker.FixedLocator([125,130,135,140,145,150,155,160])
    gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15])
    gl.xformatter   = LONGITUDE_FORMATTER
    gl.yformatter   = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':almost_black}#,'rotation': 90}
    gl.ylabel_style = {'size':10, 'color':almost_black}
    gl.ylabels_left   = True
    gl.xlabels_bottom = True

    clevs2  = [0.,0.5,1,1.5,2.,2.5,3,3.5,4,4.5,5.]
    cmap    = plt.cm.Greens#BrBG #levels=clevs2,
    print(lai_jan)
    plot2   = ax2.contourf(lons, lats, lai_jan, levels = clevs2, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
    ax2.text(0.02, 0.1, "(b)", transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=props)
    cbar = plt.colorbar(plot2, ax=ax2, ticklocation="right",orientation="horizontal", pad=0.08, aspect=25, shrink=0.40) # cax=cax, anchor=(0.0, 0.1),
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label("LAI (m$\mathregular{^{2}}$ m$\mathregular{^{-2}}$)", loc='center',size=12)
    ax2.add_patch(Polygon([[148., -36.5], [149., -36.5], [149., -35.5], [148., -35.5]],
                 closed=True,color=almost_black, fill=False))
    line = ax2.plot(x_values, y_values,c= almost_black, alpha = 1, linestyle="--", transform=ccrs.PlateCarree())
    # ================ Plotting WTD ================
    clevs1  = [1,2,3,4,5,6,7,8,9,10,11]
    cmap    = plt.cm.YlGnBu_r#BrBG_r

    # cmap    =[ (0.96862745, 0.98823529, 0.94117647, 1.        ),
    #            (0.87843137, 0.95294118, 0.85882353, 1.        ),
    #            (0.8       , 0.92156863, 0.77254902, 1.        ),
    #            (0.65882353, 0.86666667, 0.70980392, 1.        ),
    #            (0.48235294, 0.8       , 0.76862745, 1.        ),
    #            (0.30588235, 0.70196078, 0.82745098, 1.        ),
    #            (0.16862745, 0.54901961, 0.74509804, 1.        ),
    #            (0.03137255, 0.40784314, 0.6745098 , 1.        ),
    #            (0.03137255, 0.25098039, 0.50588235, 1.        ),
    #            (0.00784314, 0.21960784, 0.34509804, 1.        ),
    #            (0.03137255, 0.11372549, 0.34509804, 1.        )]

    # print("np.shape(cmap) ",np.shape(cmap))

    # 2009
    # ax3 = plt.axes(projection=ccrs.PlateCarree())
    ax3.set_extent([130,155,-44,-20])
    ax3.coastlines(resolution="50m",linewidth=1)
    ax3.add_feature(states, linewidth=.5, edgecolor="black")
    ax3.add_feature(cartopy.feature.OCEAN, edgecolor='lightgray', facecolor='lightgray')

    # Add gridlines
    gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
    gl.xlabels_top  = False
    gl.ylabels_right= False
    gl.xlines       = False
    gl.ylines       = False
    gl.xlocator     = mticker.FixedLocator([125,130,135,140,145,150,155,160])
    gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15])
    gl.xformatter   = LONGITUDE_FORMATTER
    gl.yformatter   = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':almost_black}#,'rotation': 90}
    gl.ylabel_style = {'size':10, 'color':almost_black}
    gl.ylabels_left   = True
    gl.xlabels_bottom = True

    # left - WTD
    plot3   = ax3.contourf(lons, lats, wtd_2009/1000., levels=clevs1, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
    ax3.text(0.02, 0.1, "(c)", transform=ax3.transAxes, fontsize=12, verticalalignment='top', bbox=props)
    # cbar = plt.colorbar(plot3, ax=ax3, ticklocation="right",orientation="horizontal", pad=0.03, aspect=20, shrink=0.95) # cax=cax, anchor=(0.0, 0.1),
    # cbar.set_label("m",  loc='center',size=12)
    # cbar.ax.tick_params(labelsize=10)
    ax3.add_patch(Polygon([[148., -36.5], [149., -36.5], [149., -35.5], [148., -35.5]],
                 closed=True,color=almost_black, fill=False))
    line = ax3.plot(x_values, y_values,c=almost_black, alpha = 1, linestyle="--", transform=ccrs.PlateCarree())



   # ================ Plotting WTD ================
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                        facecolor="none",
                                        name="admin_1_states_provinces_shp")

    # ax4 = plt.axes(projection=ccrs.PlateCarree())
    ax4.set_extent([130,155,-44,-20])
    ax4.coastlines(resolution="50m",linewidth=1)
    ax4.add_feature(states, linewidth=.5, edgecolor="black")
    ax4.add_feature(cartopy.feature.OCEAN, edgecolor='lightgray', facecolor='lightgray')

    # Add gridlines
    gl = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
    gl.xlabels_top  = False
    gl.ylabels_right= False
    gl.xlines       = False
    gl.ylines       = False
    gl.xlocator     = mticker.FixedLocator([125,130,135,140,145,150,155,160])
    gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15])
    gl.xformatter   = LONGITUDE_FORMATTER
    gl.yformatter   = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':almost_black}#,'rotation': 90}
    gl.ylabel_style = {'size':10, 'color':almost_black}
    gl.ylabels_left   = False
    gl.xlabels_bottom = True

    plot4   = ax4.contourf(lons, lats, wtd_2013/1000., levels=clevs1, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
    ax4.text(0.02, 0.1, "(d)", transform=ax4.transAxes, fontsize=12, verticalalignment='top', bbox=props)
    # cbar = plt.colorbar(plot4, ax=ax4, ticklocation="right",orientation="horizontal", pad=0.03, aspect=20, shrink=0.95) # cax=cax, anchor=(0.0, 0.1),
    # cbar.set_label("m",  loc='center',size=12)
    # cbar.ax.tick_params(labelsize=10)
    ax4.add_patch(Polygon([[148., -36.5], [149., -36.5], [149., -35.5], [148., -35.5]],
                 closed=True,color=almost_black, fill=False))
    line = ax4.plot(x_values, y_values,c=almost_black, alpha = 1, linestyle="--", transform=ccrs.PlateCarree())


    # 2019
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                        facecolor="none",
                                        name="admin_1_states_provinces_shp")

    # ax5 = plt.axes(projection=ccrs.PlateCarree())
    ax5.set_extent([130,155,-44,-20])
    ax5.coastlines(resolution="50m",linewidth=1)
    ax5.add_feature(states, linewidth=.5, edgecolor="black")
    ax5.add_feature(cartopy.feature.OCEAN, edgecolor='lightgray', facecolor='lightgray')

    # Add gridlines
    gl = ax5.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
    gl.xlabels_top  = False
    gl.ylabels_right= False
    gl.xlines       = False
    gl.ylines       = False
    gl.xlocator     = mticker.FixedLocator([125,130,135,140,145,150,155,160])
    gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15])
    gl.xformatter   = LONGITUDE_FORMATTER
    gl.yformatter   = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':almost_black}#,'rotation': 90}
    gl.ylabel_style = {'size':10, 'color':almost_black}
    gl.ylabels_left   = False
    gl.xlabels_bottom = True

    plot5   = ax5.contourf(lons, lats, wtd_2019/1000., levels=clevs1, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
    ax5.text(0.02, 0.1, "(e)", transform=ax5.transAxes, fontsize=12, verticalalignment='top', bbox=props)
    cbar = plt.colorbar(plot5, ax=[ax3,ax4,ax5], ticklocation="right",orientation="horizontal", pad=0.08, aspect=60, shrink=0.6) # cax=cax, anchor=(0.0, 0.1),
    cbar.set_label("WTD (m)",  loc='center',size=12)
    cbar.ax.tick_params(labelsize=10)
    ax5.add_patch(Polygon([[148., -36.5], [149., -36.5], [149., -35.5], [148., -35.5]],
                 closed=True,color=almost_black, fill=False))
    line = ax5.plot(x_values, y_values,c=almost_black, alpha = 1, linestyle="--", transform=ccrs.PlateCarree())


    fig.savefig('./plots/Fig1_spatial_map_wrf_domain_WTD_LAI', bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    #######################################################
    # Decks to run:
    #    heat_advection
    #######################################################

    file_paths  = ["/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg/LIS.CABLE.20090122-20090213_gw.nc",
                   "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2013_3Nov/ensemble_avg/LIS.CABLE.20121229-20130122_gw.nc",
                   "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2019_3Nov/ensemble_avg/LIS.CABLE.20190108-20190130_gw.nc",
                   "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2019_3Nov/gw_rst_20190103/bdy_data/lis_input.d01.nc"]
    wrf_path    =  "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg/wrfout_20090122-20090213_gw"

    spatial_map_wrf_domain(file_paths, wrf_path)
