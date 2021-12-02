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

def spatial_map_wrf_domain(file_path, wrf_path):

    # Open the NetCDF file
    encoding = 'utf-8' # Times in WRF output is btype, convert to string

    file  = Dataset(file_path)
    wrf   = Dataset(wrf_path)
    
    iveg  = file.variables['Landcover_inst'][0,:,:]
    print(np.unique(iveg))
    
    # to get lat and lon
    p     = getvar(wrf, "pressure")

    # Get the lat/lon coordinates
    lats, lons = latlon_coords(p)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(p)

    # Create the figure
    fig = plt.figure(figsize=(8,7))
    fig.subplots_adjust(hspace=0.0)
    fig.subplots_adjust(wspace=0.0)

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

    # Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([110,165,-45,-10])
        
    ax.coastlines(resolution="50m",linewidth=1)
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                         facecolor='lightgray',
                                         name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, zorder=0, edgecolor='lightgray', facecolor='lightgray') # '#fff9e6'
    ax.add_feature(cartopy.feature.LAND, edgecolor='lightgray', facecolor='lightgray')
    ax.add_feature(cartopy.feature.OCEAN, edgecolor='lightgray', facecolor='lightgray')
    ax.add_feature(cartopy.feature.COASTLINE, lw=0.5)

    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color=almost_black, linestyle='--')
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


    # ------- Plot Advection Colorfill -------
    # iveg in SE AUS : [2.0 5.0 6.0 9.0 14.0 ]
    #                  ["BEF","shrub","grass","crop","barren"]
    iveg = np.where(iveg == 2, 1, iveg)
    iveg = np.where(iveg == 5, 2, iveg)
    iveg = np.where(iveg == 6, 3, iveg)
    iveg = np.where(iveg == 9, 4, iveg)
    iveg = np.where(iveg == 14,5, iveg)
    iveg = np.where(iveg > 5  ,0, iveg)
    iveg = np.where(iveg < 1  ,0, iveg)
    
    
    cint = [-1,0.1,1.1,2.1,3.1,4.1,5.1]
    almost_black_rgb = np.array([38./255, 38./255, 38./255, 1])

    cmap = plt.cm.Paired(np.arange(12))
            # ["BEF","shrub","grass","crop","barren"]
    cmap = [cmap[0],cmap[3],cmap[2],cmap[6],cmap[7],cmap[11]] #
    cf   = ax.pcolormesh(to_np(lons), to_np(lats), iveg, levels=cint, colors=cmap,
                       transform=ccrs.PlateCarree(),rasterized=True) 

    cbar = plt.colorbar(cf, ax=ax, orientation="horizontal", pad=.07,  aspect=30, extend='neither', shrink=0.92)
    cbar.set_ticks([])
    
    ax.text(0.065, -0.15, "Water", transform=ax.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    ax.text(0.245, -0.15, "BEF", transform=ax.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    ax.text(0.39, -0.15, "Shrub", transform=ax.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    ax.text(0.55, -0.15, "Grass", transform=ax.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    ax.text(0.705, -0.15, "Crop", transform=ax.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    ax.text(0.83, -0.15, "Barren land", transform=ax.transAxes, verticalalignment='top', bbox=props, fontsize=10)
    
 
    #boundaries=[0,1,2,3,4,5], values=["Water","BEF","shrub","grass","crop","barren"])#, ticks=cint)

    fig.savefig('./plots/figures/spatial_map_wrf_domain', bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    #######################################################
    # Decks to run:
    #    heat_advection
    #######################################################

    hw_name     = "hw2019_3Nov"
    start_date  = "20190108"
    end_date    = "20190130"
    
    path        = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+hw_name+'/ensemble_avg'
    file_path   = path + '/LIS.CABLE.'+start_date+"-"+end_date+'_gw.nc'
    wrf_path    = path + '/wrfout_'+start_date+"-"+end_date+'_gw'
 
    spatial_map_wrf_domain(file_path, wrf_path)


