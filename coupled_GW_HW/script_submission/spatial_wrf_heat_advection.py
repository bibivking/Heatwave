#!/usr/bin/python

import sys
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

# import pint

def cdiff(scalar,axis=0):

    '''
    Source: https://github.com/davidmnielsen/mygrads/tree/master/mygrads
    Author: David Nielsen
    Performs the same as GrADS function cdiff()
    http://cola.gmu.edu/grads/gadoc/gradfunccdiff.html
    The scalar quantity must by 2D.
    The finite differences calculation ignores the borders, where np.nan is returned.
    '''

    # Check if 2D
    dimScalar=np.size(np.shape(scalar))
    if dimScalar != 2:
        print("Pystuff Error: scalar must have only 2 dimensions, but it has %d." %dimScalar)
        return

    # Length of each dimension
    lendim0=np.shape(scalar)[0]
    lendim1=np.shape(scalar)[1]

    # Initialize output var
    out=np.zeros(np.shape(scalar))
    out.fill(np.nan)

    # Centered finite differences
    for x in np.arange(1,lendim0-1):
        for y in np.arange(1,lendim1-1):
            if axis==0:
                out[x,y]=scalar[x+1,y]-scalar[x-1,y]
            elif axis==1:
                out[x,y]=scalar[x,y+1]-scalar[x,y-1]
            else:
                print("Pystuff Error: Invalid axis option. Must be either 0 or 1.")
                return
    return out

def cal_heat_advection(var,u,v,lat,lon):

    '''
    Calculates the advection of a scalar variable.
    lat and lon are 1D arrays.
    http://cola.gmu.edu/grads/gadoc/gradfunccdiff.html
    '''
    # latv, lonv = np.meshgrid(lat, lon, indexing='ij')

    dtr = np.pi/180
    r   = 6.371*(10**6)
    dtx = cdiff (var, axis=1)
    dty = cdiff (var, axis=0)
    dx  = cdiff(lon, axis=1)*dtr
    dy  = cdiff(lat, axis=0)*dtr
    adv = -( (u*dtx)/(np.cos(lat*dtr)*dx) + v*dty/dy )/r

    return adv

def heat_advection_metpy_rubbish(file_paths, var_name, height, time_s, time_e, var_unit=None, loc_lat=None, loc_lon=None, message=None):

    import metpy.calc as mpcalc
    from metpy.units import units

    # Open the NetCDF file
    encoding = 'utf-8' # Times in WRF output is btype, convert to string

    ncfile1  = Dataset(file_paths[0])
    ntime    = len(ncfile1.variables['Times'][:,0])
    lat      = ncfile1.variables['XLAT'][0,:,:]
    lon      = ncfile1.variables['XLONG'][0,:,:]

    time_tmp = []

    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(ncfile1.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(time_temp - datetime(2000,1,1))

    time = np.array(time_tmp)

    # to get lat and lon
    p1   = getvar(ncfile1, "pressure", timeidx=ALL_TIMES)

    # Extract the pressure, geopotential height, and wind variables
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

    # # Smooth heights and advection a little
    # # Be sure to only put in a 2D lat/lon or Y/X array for smoothing
    # Z_850 = ndimage.gaussian_filter(hght_850, sigma=3, order=0) * units.meter
    # adv = ndimage.gaussian_filter(adv, sigma=3, order=0) * units('K/sec')
    # !!!!!!!!!!!!!!

    Var1  = read_wrf_hgt_var(file_paths[0], var_name, var_unit, height, loc_lat, loc_lon)
    # Z1    = read_wrf_hgt_var(file_paths[0], "z", "m", height, loc_lat, loc_lon)
    Ua1   = read_wrf_hgt_var(file_paths[0], "ua", "m s-1", height, loc_lat, loc_lon)
    Va1   = read_wrf_hgt_var(file_paths[0], "va", "m s-1", height, loc_lat, loc_lon)

    if var_name in ['temp']:
        var1  = spital_var(time,Var1,time_s,time_e) * units.K
    else:
        scale = get_scale(var_name)
        var1  = spital_var(time,Var1,time_s,time_e)*scale

    # z1   = spital_var(time,Z1,time_s,time_e)
    ua1  = spital_var(time,Ua1,time_s,time_e) * units('m/sec')
    va1  = spital_var(time,Va1,time_s,time_e) * units('m/sec')

    # dx   = 25000
    # dy   = 25000

    print('var1 ', np.shape(var1))
    print('ua1 ', np.shape(ua1))
    print('va1 ', np.shape(va1))

    adv1 = mpcalc.advection(var1*units.K, [ua1, va1], (dx[1:,:],dy[:,1:]), x_dim=-1, y_dim=-2) * units('K/sec')
    print(np.shape(adv1))

    adv       = adv1

    # if len(file_paths) > 1:
    #     Var2  = read_wrf_hgt_var(file_paths[1], var_name, var_unit, height, loc_lat, loc_lon)
    #     Z2    = read_wrf_hgt_var(file_paths[1], "z", "m", height, loc_lat, loc_lon)
    #     Ua2   = read_wrf_hgt_var(file_paths[1], "ua", "m s-1", height, loc_lat, loc_lon)
    #     Va2   = read_wrf_hgt_var(file_paths[1], "va", "m s-1", height, loc_lat, loc_lon)
    #
    #     if var_name in ['temp']:
    #         var2  = spital_var(time,Var2,time_s,time_e)
    #     else:
    #         scale = get_scale(var_name)
    #         var2  = spital_var(time,Var2,time_s,time_e)*scale
    #
    #     z2   = spital_var(time,Z2,time_s,time_e)
    #     ua2  = spital_var(time,Ua2,time_s,time_e)
    #     va2  = spital_var(time,Va2,time_s,time_e)
    #
    #     # Calculate temperature advection using metpy function
    #     adv2 = mpcalc.advection(var2[0], [ua2[0], va2[0]], dx=25000, dy=25000, x_dim=0,y_dim=1)
    #
    #     # Calculate difference
    #     var = var2 - var1
    #     u   = ua2 - ua1
    #     v   = va2 - va1
    #     z   = z2 - z1
    #     adv = adv2- adv1
    #
    # else:
    #     # Calculate difference
    #     var = var1
    #     u   = ua1
    #     v   = va1
    #     z   = z1
    #     adv = adv1

    # Get the lat/lon coordinates
    lats, lons = latlon_coords(p1)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(p1)

    # Create the figure
    fig = plt.figure(figsize=(12,9))
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
    ax = plt.axes(projection=cart_proj)

    # Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                         facecolor="none",
                                         name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)

    # start plotting
    if loc_lat == None:
        ax.set_extent([135,155,-40,-25])
    else:
        ax.set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(p1))
    ax.set_ylim(cartopy_ylim(p1))

    # ------- Plot Advection Colorfill -------
    cint = np.arange(-8, 9)
    cf = ax.contourf(to_np(lons), to_np(lats), adv*3600.*24., cint[cint != 0],
                    extend='both', cmap='bwr', transform=ccrs.PlateCarree())
    plt.colorbar(cf, ax=ax, orientation="horizontal", pad=.05, extendrect=True, ticks=cint)
    # cb = plt.colorbar(cf, cax=cax, orientation='horizontal', extendrect=True, ticks=cint)
    # cb.set_label(r'$^{o}C/d$', size='large')


    # ------- Plot Height Contours -------
    # if len(file_paths) > 1:
    #     levels = np.arange(-2., 2., 0.1)
    # else:
    #     levels = np.arange(1500., 1600., 10.)
    #
    # contours = plt.contour(to_np(lons), to_np(lats), to_np(gaussian_filter(z,sigma=3)),
    #                        levels = levels, colors="black", linewidths=1.5, linestyles='solid',
    #                        transform=ccrs.PlateCarree())
    # plt.clabel(contours, inline=1, inline_spacing=10, fontsize=6, fmt="%d",rightside_up=True, use_clabeltext=True)

    # ------- Plot Var Contours ----------
    # if len(file_paths) > 1:
    #     max_val = np.abs(np.nanmax(var))
    #     min_val = np.abs(np.nanmin(var))
    #     max_range = np.maximum(max_val,min_val)
    #     print(max_val)
    #     print(min_val)
    #     print(max_range)
    #     levels = np.linspace(max_range*(-1.),max_range,num=20)
    # else:
    #     levels = np.arange(np.nanmin(var), np.nanmax(var), 20)
    #
    # var_contours = plt.contour(to_np(lons), to_np(lats), to_np(var),levels = levels,
    #                            colors='grey', linewidths=1.25, linestyles='dashed',
    #                            transform=ccrs.PlateCarree()) #,"jet" #“rainbow”#"coolwarm" , cmap=get_cmap("bwr"),extend='both'
    # # plt.colorbar(var_contours, ax=ax, orientation="horizontal", pad=.05)
    # plt.clabel(var_contours, fontsize=6, inline=1, inline_spacing=10, fmt='%d',
    #             rightside_up=True, use_clabeltext=True)

    # -------- Plot wind barbs --------
    # if len(file_paths) > 1:
    #     ax.quiver(to_np(lons[::3,::3]), to_np(lats[::3,::3]), to_np(u[::3, ::3]),to_np(v[::3, ::3]),
    #          scale=20., pivot='middle', transform=ccrs.PlateCarree()) # width=0.0002,
    # else:
    #     ax.quiver(to_np(lons[::3,::3]), to_np(lats[::3,::3]), to_np(u[::3, ::3]),to_np(v[::3, ::3]),
    #          scale=100., pivot='middle', transform=ccrs.PlateCarree()) # width=0.0002,


    if var_unit == None:
        plt.title(str(height)+"hPa, "+var_name)
    else:
        plt.title(str(height)+"hPa, "+var_name+" advection ("+var_unit+"), Geopotential Height (gpm), and Barbs (m s-1)")

    if message == None:
        message = var_name+'_'+str(height)+"hPa"
    else:
        message = message+"_"+var_name+'_'+str(height)+"hPa"

    fig.savefig('./plots/5Nov/heat_advection/3Nov/spatial_map_wrf_heat_advection_'+message , bbox_inches='tight', pad_inches=0.1)

def heat_advection(file_paths, var_name, height, time_s, time_e, seconds=None, var_unit=None, loc_lat=None, loc_lon=None, message=None):

    # Open the NetCDF file
    encoding = 'utf-8' # Times in WRF output is btype, convert to string

    ncfile1  = Dataset(file_paths[0])
    ntime    = len(ncfile1.variables['Times'][:,0])
    lat      = ncfile1.variables['XLAT'][0,:,:]
    lon      = ncfile1.variables['XLONG'][0,:,:]

    time_tmp = []

    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(ncfile1.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1))

    time = np.array(time_tmp)

    # to get lat and lon
    p1   = getvar(ncfile1, "pressure", timeidx=ALL_TIMES)

    Var1  = read_wrf_hgt_var(file_paths[0], var_name, var_unit, height, loc_lat, loc_lon)
    Ua1   = read_wrf_hgt_var(file_paths[0], "ua", "m s-1", height, loc_lat, loc_lon)
    Va1   = read_wrf_hgt_var(file_paths[0], "va", "m s-1", height, loc_lat, loc_lon)

    if var_name in ['temp']:
        var1  = spital_var(time,Var1,time_s,time_e,seconds)
    else:
        scale = get_scale(var_name)
        var1  = spital_var(time,Var1,time_s,time_e,seconds)*scale

    # z1   = spital_var(time,Z1,time_s,time_e)
    ua1  = spital_var(time,Ua1,time_s,time_e,seconds)
    va1  = spital_var(time,Va1,time_s,time_e,seconds)

    adv1 = cal_heat_advection(var1,ua1,va1,lat,lon)

    print(np.shape(adv1))

    if len(file_paths) > 1:
        Var2  = read_wrf_hgt_var(file_paths[1], var_name, var_unit, height, loc_lat, loc_lon)
        Ua2   = read_wrf_hgt_var(file_paths[1], "ua", "m s-1", height, loc_lat, loc_lon)
        Va2   = read_wrf_hgt_var(file_paths[1], "va", "m s-1", height, loc_lat, loc_lon)

        if var_name in ['temp']:
            var2  = spital_var(time,Var2,time_s,time_e,seconds)
        else:
            scale = get_scale(var_name)
            var2  = spital_var(time,Var2,time_s,time_e,seconds)*scale

        # z2   = spital_var(time,Z2,time_s,time_e)
        ua2  = spital_var(time,Ua2,time_s,time_e,seconds)
        va2  = spital_var(time,Va2,time_s,time_e,seconds)

        adv2 = cal_heat_advection(var2,ua2,va2,lat,lon)

        # Calculate difference
        adv = (adv2 - adv1)*3600.*24.
        ua  = ua2  - ua1
        va  = va2  - va1
    else:
        adv = adv1*3600.*24.
        ua  = ua1
        va  = va1

    # Get the lat/lon coordinates
    lats, lons = latlon_coords(p1)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(p1)

    # Create the figure
    fig = plt.figure(figsize=(12,9))
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
    ax.set_extent([125,160,-45,-20])
    ax.coastlines(resolution="50m",linewidth=1)
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                         facecolor="none",
                                         name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")

    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
    gl.xlabels_top  = False
    gl.ylabels_right= False
    gl.xlines       = True
    gl.xlocator     = mticker.FixedLocator([125,130,135,140,145,150,155,160])
    gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20])
    gl.xformatter   = LONGITUDE_FORMATTER
    gl.yformatter   = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':'black'}#,'rotation': 90}
    gl.ylabel_style = {'size':10, 'color':'black'}

    # Set the map bounds
    # ax.set_xlim(cartopy_xlim(p1))
    # ax.set_ylim(cartopy_ylim(p1))

    # # start plotting
    # if loc_lat == None:
    #     ax.set_extent([135,155,-40,-25])
    # else:
    #     ax.set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])

    # ------- Plot Advection Colorfill -------

    if len(file_paths) > 1:
        cint = [-6,-5,-4,-3,-2,-1,1,2,3,4,5,6]
    else:
        cint = [-30,-25,-20,-15,-10,-5,5,10,15,20,25,30]

    cf = ax.contourf(to_np(lons), to_np(lats), adv, levels=cint,
                    extend='both', cmap='bwr', transform=ccrs.PlateCarree())
                    #
    plt.colorbar(cf, ax=ax, orientation="horizontal", pad=.05, extendrect=True, ticks=cint)
    # cb.set_label(r'$^{o}C/d$', size='large')


    # ------- Plot Height Contours -------
    # if len(file_paths) > 1:
    #     levels = np.arange(-2., 2., 0.1)
    # else:
    #     levels = np.arange(1500., 1600., 10.)
    #
    # contours = plt.contour(to_np(lons), to_np(lats), to_np(gaussian_filter(z,sigma=3)),
    #                        levels = levels, colors="black", linewidths=1.5, linestyles='solid',
    #                        transform=ccrs.PlateCarree())
    # plt.clabel(contours, inline=1, inline_spacing=10, fontsize=6, fmt="%d",rightside_up=True, use_clabeltext=True)

    # -------- Plot wind barbs --------
    if len(file_paths) > 1:
        scale = 5.
    else:
        scale = 10.

    q = ax.quiver(to_np(lons[::3,::3]), to_np(lats[::3,::3]), to_np(ua[::3, ::3]),to_np(va[::3, ::3]),
            angles='xy', scale_units='xy' , scale=scale, pivot='middle', transform=ccrs.PlateCarree()) # width=0.0002,
    ax.quiverkey(q, X=0.1, Y=0.9, U=scale, label=str(scale)+' m/s', labelpos='E')

    if var_unit == None:
        plt.title(str(height)+"hPa, "+var_name)
    else:
        plt.title(str(height)+"hPa, "+var_name+" advection (K d-1) and Barbs (m s-1)")

    if message == None:
        message = var_name+'_'+str(height)+"hPa"
    else:
        message = message+"_"+var_name+'_'+str(height)+"hPa"

    fig.savefig('./plots/5Nov/heat_advection/3Nov/spatial_map_wrf_heat_advection_'+message , bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":

    #######################################################
    # Decks to run:
    #    heat_advection
    #######################################################

    hw_name           = "hw2009_3Nov"
    var_name          = 'temp'
    var_unit          = 'K'#'degC'
    height            = 850
    seconds           = [6.*60.*60.,18.*60.*60.]

    # 2009
    if hw_name == "hw2009_3Nov":
        start_date= "20090122"
        end_date  = "20090213"
        rst_dates = ["20090117","20090118","20090119","20090120","20090121" ]
        time_s = datetime(2009,1,28,0,0,0,0)
        time_e = datetime(2009,2,8,23,59,0,0)
        # Time_s = datetime(2009,1,22,0,0,0,0)
        # Time_e = datetime(2009,2,13,23,59,0,0)

    elif hw_name == "hw2013_3Nov":
        start_date= "20121229"
        end_date  = "20130122"
        rst_dates = ["20121224","20121225","20121226","20121227","20121228" ]
        time_s = datetime(2013,1,4,0,0,0,0)
        time_e = datetime(2013,1,18,23,59,0,0)
        # Time_s = datetime(2012,12,29,0,0,0,0)
        # Time_e = datetime(2013,1,22,23,59,0,0)

    elif hw_name == "hw2019_3Nov":
        start_date= "20190108"
        end_date  = "20190130"
        rst_dates = ["20190103","20190104","20190105","20190106","20190107" ]

        time_s = datetime(2019,1,14,0,0,0)
        time_e = datetime(2019,1,26,23,59,0,0)
        # Time_s = datetime(2019,1,8,14,0,0,0)
        # Time_e = datetime(2019,1,30,0,0,0,0)


    cpl_atmo_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+hw_name+'/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_'+start_date+"-"+end_date+'_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_'+start_date+"-"+end_date+'_fd'  # atmo output of wrf-cable run

    file_paths        = [cpl_atmo_file_fd,cpl_atmo_file_gw] #cpl_atmo_file_fd, cpl_atmo_file_fd,cpl_atmo_file_fd, cpl_atmo_file_gw


    # tot_day   = (Time_e-Time_s).days + 1

    # for i in np.arange(tot_day):
        # time_s = Time_s + timedelta(days=int(i))
        # time_e = Time_s + timedelta(days=int(i+1)) - timedelta(seconds=1)
    if len(file_paths) > 1:
        message = "Couple_GW-FD_"+str(time_s)+"_"+str(time_e)
    else:
        message = "Couple_GW_"+str(time_s)+"_"+str(time_e)

    heat_advection(file_paths, var_name, height, time_s, time_e, seconds, var_unit, message=message)
