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
from cartopy.feature import NaturalEarthFeature
import metpy.calc as mpcalc
from metpy.units import units
import scipy.ndimage as ndimage
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                        cartopy_ylim, latlon_coords, ALL_TIMES)
from common_utils import *

def heat_advection(file_paths, var_name, height, time_s, time_e, var_unit=None, loc_lat=None, loc_lon=None, message=None):

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
    p1       = getvar(ncfile1, "pressure", timeidx=ALL_TIMES)

    # Extract the pressure, geopotential height, and wind variables
    # dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)
    # print(dx)
    # # Smooth heights and advection a little
    # # Be sure to only put in a 2D lat/lon or Y/X array for smoothing
    # Z_850 = ndimage.gaussian_filter(hght_850, sigma=3, order=0) * units.meter
    # adv = ndimage.gaussian_filter(adv, sigma=3, order=0) * units('K/sec')
    # !!!!!!!!!!!!!!

    Var1  = read_wrf_hgt_var(file_paths[0], var_name, var_unit, height, loc_lat, loc_lon)
    Z1    = read_wrf_hgt_var(file_paths[0], "z", "m", height, loc_lat, loc_lon)
    Ua1   = read_wrf_hgt_var(file_paths[0], "ua", "m s-1", height, loc_lat, loc_lon)
    Va1   = read_wrf_hgt_var(file_paths[0], "va", "m s-1", height, loc_lat, loc_lon)

    if var_name in ['temp']:
        var1  = spital_var(time,Var1,time_s,time_e)
    else:
        scale = get_scale(var_name)
        var1  = spital_var(time,Var1,time_s,time_e)*scale
    
    z1   = spital_var(time,Z1,time_s,time_e)
    ua1  = spital_var(time,Ua1,time_s,time_e)
    va1  = spital_var(time,Va1,time_s,time_e)

    # Calculate temperature advection using metpy function
    # dx   = 25000
    # dy   = 25000
    adv1 = mpcalc.advection(var1, [ua1, va1], dx=25000, dy=25000, x_dim=-2,y_dim=-1) #, dim_order='yx')  #  * units.kelvin * units('K/sec')
                        
    if len(file_paths) > 1:
        Var2  = read_wrf_hgt_var(file_paths[1], var_name, var_unit, height, loc_lat, loc_lon)
        Z2    = read_wrf_hgt_var(file_paths[1], "z", "m", height, loc_lat, loc_lon)
        Ua2   = read_wrf_hgt_var(file_paths[1], "ua", "m s-1", height, loc_lat, loc_lon)
        Va2   = read_wrf_hgt_var(file_paths[1], "va", "m s-1", height, loc_lat, loc_lon)

        if var_name in ['temp']:
            var2  = spital_var(time,Var2,time_s,time_e)
        else:
            scale = get_scale(var_name)
            var2  = spital_var(time,Var2,time_s,time_e)*scale
        
        z2   = spital_var(time,Z2,time_s,time_e)
        ua2  = spital_var(time,Ua2,time_s,time_e)
        va2  = spital_var(time,Va2,time_s,time_e)

        # Calculate temperature advection using metpy function
        adv2 = mpcalc.advection(var2, [ua2, va2], dx=25000, dy=25000, x_dim=-2,y_dim=-1) 

        # Calculate difference
        var = var2 - var1
        u   = ua2 - ua1
        v   = va2 - va1
        z   = z2 - z1
        adv = adv2- adv1
        
    else:
        # Calculate difference
        var = var1
        u   = ua1
        v   = va1
        z   = z1
        adv = adv1

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
    cb.set_label(r'$^{o}C/d$', size='large')


    # ------- Plot Height Contours -------
    if len(file_paths) > 1:
        levels = np.arange(-2., 2., 0.1)
    else:
        levels = np.arange(1500., 1600., 10.)
    
    contours = plt.contour(to_np(lons), to_np(lats), to_np(gaussian_filter(z,sigma=3)), 
                           levels = levels, colors="black", linewidths=1.5, linestyles='solid',
                           transform=ccrs.PlateCarree())
    plt.clabel(contours, inline=1, inline_spacing=10, fontsize=6, fmt="%d",rightside_up=True, use_clabeltext=True)

    # ------- Plot Var Contours ----------
    if len(file_paths) > 1:
        max_val = np.abs(np.nanmax(var))
        min_val = np.abs(np.nanmin(var))
        max_range = np.maximum(max_val,min_val)
        print(max_val)
        print(min_val)
        print(max_range)
        levels = np.linspace(max_range*(-1.),max_range,num=20)
    else:
        levels = np.arange(np.nanmin(var), np.nanmax(var), 20)

    var_contours = plt.contour(to_np(lons), to_np(lats), to_np(var),levels = levels, 
                               colors='grey', linewidths=1.25, linestyles='dashed',
                               transform=ccrs.PlateCarree()) #,"jet" #“rainbow”#"coolwarm" , cmap=get_cmap("bwr"),extend='both'
    # plt.colorbar(var_contours, ax=ax, orientation="horizontal", pad=.05)
    plt.clabel(var_contours, fontsize=6, inline=1, inline_spacing=10, fmt='%d',
                rightside_up=True, use_clabeltext=True)

    # -------- Plot wind barbs -------- 
    if len(file_paths) > 1:
        ax.quiver(to_np(lons[::3,::3]), to_np(lats[::3,::3]), to_np(u[::3, ::3]),to_np(v[::3, ::3]),
             scale=20., pivot='middle', transform=ccrs.PlateCarree()) # width=0.0002,
    else:
        ax.quiver(to_np(lons[::3,::3]), to_np(lats[::3,::3]), to_np(u[::3, ::3]),to_np(v[::3, ::3]),
             scale=100., pivot='middle', transform=ccrs.PlateCarree()) # width=0.0002,


    if var_unit == None:
        plt.title(str(height)+"hPa, "+var_name)
    else:
        plt.title(str(height)+"hPa, "+var_name+" advection ("+var_unit+"), Geopotential Height (gpm), and Barbs (m s-1)")

    if message == None:
        message = var_name+'_'+str(height)+"hPa"
    else:
        message = message+"_"+var_name+'_'+str(height)+"hPa"

    fig.savefig('./plots/5Nov/hw_anaylsis/3Nov/spatial_map_wrf_heat_advection_'+message , bbox_inches='tight', pad_inches=0.1)

    # # Helper function for finding proper time variable
    # def find_time_var(var, time_basename='time'):
    #     for coord_name in var.coordinates.split():
    #         if coord_name.startswith(time_basename):
    #             return coord_name
    #     raise ValueError('No time variable found for ' + var.name)

    # # Create NCSS object to access the NetcdfSubset
    # # Data from NCEI GFS 0.5 deg Analysis Archive

    # base_url = 'https://www.ncei.noaa.gov/thredds/ncss/grid/gfs-g4-anl-files/'
    # dt = datetime(2017, 4, 5, 12)
    # ncss = NCSS('{}{dt:%Y%m}/{dt:%Y%m%d}/gfsanl_4_{dt:%Y%m%d}_'
    #             '{dt:%H}00_000.grb2'.format(base_url, dt=dt))

    # # Create lat/lon box for location you want to get data for
    # query = ncss.query().time(dt)
    # query.lonlat_box(north=65, south=15, east=310, west=220)
    # query.accept('netcdf')

    # # Request data for vorticity
    # query.variables('Geopotential_height_isobaric', 'Temperature_isobaric',
    #                 'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric')
    # data = ncss.get_data(query)

    # # Pull out variables you want to use
    # hght_var = data.variables['Geopotential_height_isobaric']
    # temp_var = data.variables['Temperature_isobaric']
    # u_wind_var = data.variables['u-component_of_wind_isobaric']
    # v_wind_var = data.variables['v-component_of_wind_isobaric']
    # time_var = data.variables[find_time_var(temp_var)]
    # lat_var = data.variables['lat']
    # lon_var = data.variables['lon']

    # # Get actual data values and remove any size 1 dimensions
    # lat = lat_var[:].squeeze()
    # lon = lon_var[:].squeeze()
    # hght = hght_var[:].squeeze()
    # temp = temp_var[:].squeeze() * units.K
    # u_wind = units('m/s') * u_wind_var[:].squeeze()
    # v_wind = units('m/s') * v_wind_var[:].squeeze()

    # # Convert number of hours since the reference time into an actual date
    # time = num2date(time_var[:].squeeze(), time_var.units)

    # lev_850 = np.where(data.variables['isobaric'][:] == 850*100)[0][0]
    # hght_850 = hght[lev_850]
    # temp_850 = temp[lev_850]
    # u_wind_850 = u_wind[lev_850]
    # v_wind_850 = v_wind[lev_850]

    # # Combine 1D latitude and longitudes into a 2D grid of locations
    # lon_2d, lat_2d = np.meshgrid(lon, lat)
    # # Gridshift for barbs
    # lon_2d[lon_2d > 180] = lon_2d[lon_2d > 180] - 360
    # # Begin data calculations
    # # Use helper function defined above to calculate distance
    # # between lat/lon grid points
    # dx, dy = mpcalc.lat_lon_grid_deltas(lon_var, lat_var)

    # # Calculate temperature advection using metpy function
    
    # adv = mpcalc.advection(temp_850 * units.kelvin, [u_wind_850, v_wind_850],
    #                     (dx, dy), dim_order='yx') * units('K/sec')

    # # Smooth heights and advection a little
    # # Be sure to only put in a 2D lat/lon or Y/X array for smoothing
    # Z_850 = ndimage.gaussian_filter(hght_850, sigma=3, order=0) * units.meter
    # adv = ndimage.gaussian_filter(adv, sigma=3, order=0) * units('K/sec')

    # # Begin map creation
    # # Set Projection of Data
    # datacrs = ccrs.PlateCarree()

    # # Set Projection of Plot
    # plotcrs = ccrs.LambertConformal(central_latitude=[30, 60], central_longitude=-100)

    # # Create new figure
    # fig = plt.figure(figsize=(11, 8.5))
    # gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
    #                     hspace=0.01, wspace=0.01)

    # # Add the map and set the extent
    # ax = plt.subplot(gs[0], projection=plotcrs)
    # plt.title('850mb Temperature Advection for {0:%d %B %Y %H:%MZ}'.format(time), fontsize=16)
    # ax.set_extent([235., 290., 20., 55.])

    # # Add state/country boundaries to plot
    # ax.add_feature(cfeature.STATES)
    # ax.add_feature(cfeature.BORDERS)

    # !!!!!!!!!!!!!
    # # Plot Height Contours
    # clev850 = np.arange(900, 3000, 30)
    # cs = ax.contour(lon_2d, lat_2d, Z_850, clev850, colors='black', linewidths=1.5,
    #                 linestyles='solid', transform=datacrs)
    # plt.clabel(cs, fontsize=10, inline=1, inline_spacing=10, fmt='%i',
    #         rightside_up=True, use_clabeltext=True)

    # # Plot Temperature Contours
    # clevtemp850 = np.arange(-20, 20, 2)
    # cs2 = ax.contour(lon_2d, lat_2d, temp_850.to(units('degC')), clevtemp850,
    #                 colors='grey', linewidths=1.25, linestyles='dashed',
    #                 transform=datacrs)
    # plt.clabel(cs2, fontsize=10, inline=1, inline_spacing=10, fmt='%i',
    #         rightside_up=True, use_clabeltext=True)

    # # Plot Colorfill of Temperature Advection
    # cint = np.arange(-8, 9)
    # cf = ax.contourf(lon_2d, lat_2d, 3*adv.to(units('delta_degC/hour')), cint[cint != 0],
    #                 extend='both', cmap='bwr', transform=datacrs)
    # cax = plt.subplot(gs[1])
    # cb = plt.colorbar(cf, cax=cax, orientation='horizontal', extendrect=True, ticks=cint)
    # cb.set_label(r'$^{o}C/3h$', size='large')

    # # Plot Wind Barbs
    # ax.barbs(lon_2d, lat_2d, u_wind_850.magnitude, v_wind_850.magnitude,
    #         length=6, regrid_shape=20, pivot='middle', transform=datacrs)

    # plt.show()
    # !!!!!!!!!!!

    # vort_adv = mpcalc.advection(vor, [u, v], (dx, dy), dim_order='yx') * 1e9  上面都是对的，就这句出错


    # def heat_content():



if __name__ == "__main__":

    #######################################################
    # Decks to run:
    #    heat_advection
    #######################################################

    cpl_atmo_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_3Nov/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_20090122-20090213_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_20090122-20090213_fd'  # atmo output of wrf-cable run

    file_paths        = [cpl_atmo_file_fd,cpl_atmo_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw

    var_name          = 'temp'
    var_unit          = 'degC'
    height            = 850

    # for i in np.arange(0,23):
    time_s = datetime(2009,1,29,14,0,0,0) #+ timedelta(days=int(i))
    time_e = datetime(2009,1,30,13,59,0,0) #+ timedelta(days=int(i))

    if len(file_paths) > 1:
        message = "Couple_GW-FD_"+str(time_s)
    else:
        message = "Couple_GW_"+str(time_s)

    heat_advection(file_paths, var_name, height, time_s, time_e, var_unit, message=message) #  loc_lat=loc_lat, loc_lon=loc_lat,
