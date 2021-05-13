from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

# def plot_Tair_Wind(ftair, fuas, fvas):

# date = ''
path = '/g/data/ub4/erai/netcdf/3hr/atmos/oper_fc_sfc/v01'
ftair = path + '/tas/tas_3hrs_ERAI_historical_fc-sfc_20130101_20130131.nc'
fuas  = path + '/uas/uas_3hrs_ERAI_historical_fc-sfc_20130101_20130131.nc'
fvas  = path + '/vas/vas_3hrs_ERAI_historical_fc-sfc_20130101_20130131.nc'

# Open the NetCDF4 file (add a directory path if necessary) for reading:
tas = Dataset(ftair, mode='r')
uas = Dataset(fuas, mode='r')
vas = Dataset(fvas, mode='r')
# Run the following cell to see the MERRA2 metadata. This line will print attribute and variable information. From the 'variables(dimensions)' list, choose which variable(s) to read in below:
# print(tas)
# Read in variables:
# longitude and latitude
lons = tas.variables['lon']
lats = tas.variables['lat']
lon, lat = np.meshgrid(lons, lats)

# 2-meter air temperatur K
T2M = tas.variables['tas']

# 10-meter eastward wind m/s
U10M = uas.variables['uas']

# 10-meter northward wind m/s
V10M = vas.variables['vas']

# # Replace _FillValues with NaNs:
# U10M_nans = U10M[:]
# V10M_nans = V10M[:]
# _FillValueU10M = U10M._FillValue
# _FillValueV10M = V10M._FillValue
# U10M_nans[U10M_nans == _FillValueU10M] = np.nan
# V10M_nans[V10M_nans == _FillValueV10M] = np.nan

# # Calculate wind speed:
# ws = np.sqrt(U10M_nans**2+V10M_nans**2)
# # Calculate wind direction in radians:
# ws_direction = np.arctan2(V10M_nans,U10M_nans)

# for i in np.arange(31):
    # NOTE: the MERRA-2 file contains hourly data for 24 hours (t=24). To get the daily mean wind speed, take the average of the hourly wind speeds:
    # T2M_daily_avg = np.nanmean(T2M[i*8:(i+1)*8,:,:], axis=0) -273.15 # i*8:(i+1)*8
    # print(T2M_daily_avg)
    # # NOTE: To calculate the average wind direction correctly it is important to use the 'vector average' as atan2(<v>,<u>) where <v> and <u> are the daily average component vectors, rather than as mean of the individual wind vector direction angle.  This avoids a situation where averaging 1 and 359 = 180 rather than the desired 0.
    # U10M_daily_avg = np.nanmean(U10M[i*8:(i+1)*8,:,:], axis=0)
    # V10M_daily_avg = np.nanmean(V10M[i*8:(i+1)*8,:,:], axis=0)

    # ws_daily_avg_direction = np.arctan2(V10M_daily_avg, U10M_daily_avg)

    #Plot Global MERRA-2 Wind Speed
    # Set the figure size, projection, and extent
    # fig = plt.figure(figsize=(8,4))
    # ax = plt.axes(projection=ccrs.Robinson())
    # ax.set_global()
    # ax.coastlines(resolution="110m",linewidth=1)
    # ax.gridlines(linestyle='--',color='black')
    # # Plot windspeed: set contour levels, then draw the filled contours and a colorbar
    # clevs = np.arange(0,19,1)
    # plt.contourf(lon, lat, ws_daily_avg, clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.jet)
    #
    # plt.title('MERRA-2 Daily Average 2-meter Wind Speed, 1 June 2010', size=14)
    # cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    # cb.set_label('m/s',size=12,rotation=0,labelpad=15)
    # cb.ax.tick_params(labelsize=10)
    # plt.savefig('Tair.png',dpi=1200)
    # plt.show()

for i in np.arange(31*8):
    # The filled contours show the wind speed. The "quiver" function is used to overlay arrows to show the wind direction. The length of the arrows is determined by the wind speed.
    # Set the figure size, projection, and extent
    fig = plt.figure(figsize=(9,5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([110,155,-45,-10])
    ax.coastlines(resolution="50m",linewidth=1)
    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = True
    gl.xlocator = mticker.FixedLocator([110,115,120,125,130,135,140,145,150,155])
    gl.ylocator = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15,-10])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':'black'}
    gl.ylabel_style = {'size':10, 'color':'black'}
    # Plot windspeed
    clevs = np.arange(0,44,2)
    plt.contourf(lon, lat, T2M[i,:,:]-273.15, clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.jet) # T2M_daily_avg
    plt.title('ERA5 2m Air Temperture ' + str(i), size=16)
    cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    cb.set_label('K',size=14,rotation=0,labelpad=15)
    cb.ax.tick_params(labelsize=10)
    # Overlay wind vectors
    qv = plt.quiver(lon, lat, U10M[i,:,:], V10M[i,:,:], scale=420, color='k')
    plt.savefig('Tair_Wind_201301-'+str(i)+'.png',dpi=1200)
