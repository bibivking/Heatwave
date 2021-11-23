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
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_range_diff
from common_utils import *


def calc_diff(var,axis=0):

    '''

    '''

    nx = np.shape(var)[0]-1 
    ny = np.shape(var)[1]-1
    nz = np.shape(var)[2]

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

def air_density(T,QV,P):

    # source: https://mailman.ucar.edu/pipermail/ncl-talk/2019-August/015431.html

    '''
    Ideal gas law:  p = (density)*RD*(temperature)
                        Hence:
                        density = p/(RD*temperature)
    p   => Pa
    qv  => kg/kg  :QVAPOR
    qc  => kg/kg  :QCLOUD
    t   => K
    '''

    RD      = 287.04         # [J/(K-kg)] Gas Constant dry air

    TV      = T*(1+0.609*QV)
    density = P/(RD*TV)

    return density

def calc_atmospheric_heat_content(T,Density,dx,dy,dz):

    """
    Quantify ocean heat content upto certain depths.
    param T: interpolation of potential temperature on the vector grid (v grid)    [z, jj, ji]
    param mask: land-sea mask on vector grid               [z, jj, ji]
    param e1t: width of the grid box in zonal direction on T grid    [jj, ji]
    param e2t: height of the grid box in zonal direction on T grid   [jj, ji]
    param e3t_0: depth of each grid box                    [z]
    param e3t_adjust: adjustment of grid box depth based on the definition of partial grid  [z, jj, ji]

    return: arrays of OHC upto certain depth
    rtype: numpy arrays 
    """

    # define the constant:
    constant = {
                'g' : 9.80616,      # gravititional acceleration [m / s2]
                'cp': 1006,         # heat capacity of air [J/(Kg*C)]
                'rho': 1027,        # sea water density [Kg/m3]
                }

    nlat = len(dx[:,0])
    nlon = len(dx[0,:])
    nelv = len(dz[:,0,0])

    # calculate heat flux at each grid point
    AHC = np.zeros((nelv,nlat,nlon),dtype=float)

    for i in np.arange(nelv):
        AHC[i,:,:] = Density[i,:,:]*constant['cp']*T[i,:,:]*dx*dy*dz[i]
        AHC_sum = np.sum(AHC,0) /1e+12
        
    return AHC_sum

def spatial_map_heat_content(file_paths, time_s, time_e, seconds=None, message=None):

    # Open the NetCDF file
    encoding = 'utf-8' # Times in WRF output is btype, convert to string

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    file1    = Dataset(file_paths[0], mode='r')
    lon      = file1.variables['XLONG'][0,:,:]
    lat      = file1.variables['XLAT'][0,:,:]

    ntime    = len(file1.variables['Times'][:,0])
    nlat     = len(lat[:,0])
    nlon     = len(lat[0,:])

    time_tmp = []

    for i in np.arange(ntime):
        time_temp = datetime.strptime(str(file1.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
        time_tmp.append(UTC_to_AEST(time_temp) - datetime(2000,1,1))

    time = np.array(time_tmp)
    
    PBLH1  = read_wrf_hgt_var(file_paths[0], "PBLH", "m")  
    Z1     = read_wrf_hgt_var(file_paths[0], "height_agl", "m") 
    T1     = read_wrf_hgt_var(file_paths[0], "temp", "K") 
    QV1    = read_wrf_hgt_var(file_paths[0], "QVAPOR") # kg kg-1
    P1     = read_wrf_hgt_var(file_paths[0], "p", "Pa")


    pblh1  = spital_var(time,PBLH1,time_s,time_e,seconds)
    z1     = spital_var(time,Z1,time_s,time_e,seconds)
    t1     = spital_var(time,T1,time_s,time_e,seconds)
    qv1    = spital_var(time,QV1,time_s,time_e,seconds)
    p1     = spital_var(time,P1,time_s,time_e,seconds)

    print(np.shape(pblh1))
    print(np.shape(z1))
    print(np.shape(t1))
    print(np.shape(qv1))
    print(np.shape(p1))

    pblh1    
    # for x in np.arange(nlat):
    #     for y in np.arange(nlon):
    #         if x == 0 or x == nlat-1:
    #             if y == 
    #         z1[:,:,:] 
    # density1 = air_density(t1,qv1,p1)
    # atmo_hc1 = calc_atmospheric_heat_content(t1,density1,dx,dy,dz)

    # if len(file_paths) > 1:
    #     file2 = Dataset(file_paths[1], mode='r')
    #     T2    = read_wrf_hgt_var(file_paths[1], "temp", "K") 
    #     QV2   = read_wrf_hgt_var(file_paths[1], "QVAPOR") # kg kg-1
    #     P2    = read_wrf_hgt_var(file_paths[1], "p", "Pa")

    #     t2    = spital_var(time,T2,time_s,time_e,seconds)
    #     qv2   = spital_var(time,QV1,time_s,time_e,seconds)
    #     p2    = spital_var(time,P1,time_s,time_e,seconds)
    #     density2 = air_density(t2,qv2,p2)
    #     atmo_hc2 = calc_atmospheric_heat_content(t2,density2,dx,dy,dz)

    #     atmo_hc  = atmo_hc2 - atmo_hc1
    # else:
    #     atmo_hc  = atmo_hc1


    # # Make plots
    # fig = plt.figure(figsize=(7,5))
    # ax = plt.axes(projection=ccrs.PlateCarree())
    # # ax.set_extent([135,155,-40,-25])
    # ax.coastlines(resolution="50m",linewidth=1)

    # # Add gridlines
    # gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    # gl.xlines = True
    # gl.xlocator     = mticker.FixedLocator([110,115,120,125,130,135,140,145,150,155])
    # gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15,-10])
    # gl.xformatter   = LONGITUDE_FORMATTER
    # gl.yformatter   = LATITUDE_FORMATTER
    # gl.xlabel_style = {'size':10, 'color':'black'}
    # gl.ylabel_style = {'size':10, 'color':'black'}

    # if get_reverse_colormap(var_name) == None:
    #     cmap = plt.cm.seismic
    # elif get_reverse_colormap(var_name) == True:
    #     cmap = plt.cm.seismic_r
    # else:
    #     cmap = plt.cm.seismic_r

    # if len(file_paths) > 1:
    #     if var_name == 'EF':
    #         clevs = np.linspace(-0.5,0.5, num=21)
    #     elif var_name == 'Qair_f_inst':
    #         clevs = np.linspace(-2,2, num=21)
    #         var = var*1000. # kg/kg -> g/kg
    #     elif var_name == 'Tair_f_inst':
    #         clevs = np.linspace(-2,2, num=21)
    #     else:
    #         clevs = np.linspace(-50.,50., num=21)
    # else:
    #     clevs = np.arange(np.nanmin(var), np.nanmax(var), 21)

    # plt.contourf(lon, lat, var, levels=clevs[clevs!=0], transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
    # cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

    # plt.title(var_name, size=16)
    # cb.ax.tick_params(labelsize=10)

    # if message == None:
    #     message = var_name
    # else:
    #     message = message+"_"+var_name

    # plt.savefig('./plots/5Nov/land_var/3Nov/spatial_map_'+message+'.png',dpi=300)
 
if __name__ == "__main__":
    
    # =============================== Operation ================================
    case_name    = "hw2009_3Nov" # "hw2013_3Nov"# "hw2019_3Nov"

    if case_name == "hw2009_3Nov":
        period = "20090122-20090213"
        time_s = datetime(2009,1,28,0,0,0,0)
        time_e = datetime(2009,2,8,23,59,0,0)

    elif case_name == "hw2013_3Nov":
        period = "20121229-20130122"
        time_s = datetime(2013,1,4,0,0,0,0)
        time_e = datetime(2013,1,18,23,59,0,0)

    elif case_name == "hw2019_3Nov":
        period = "20190108-20190130"
        time_s = datetime(2019,1,14,0,0,0)
        time_e = datetime(2019,1,26,23,59,0,0)

    cpl_atmo_file     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_'+period+'_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_'+period+'_fd'  # atmo output of wrf-cable run

    file_paths        = [cpl_atmo_file_fd, cpl_atmo_file_gw] # cpl_atmo_file_fd, cpl_atmo_file_gw
    seconds           = None #[18.*60.*60.,6.*60.*60.]

    if len(file_paths) > 1:
        message = "heat_content_GW-FD_"+str(time_s)+"-"+str(time_e)
    else:
        message = "heat_content_GW_"+str(time_s)+"-"+str(time_e)

    spatial_map_heat_content(file_paths, time_s, time_e, seconds=seconds, message=message)
