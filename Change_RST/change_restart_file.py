#!/usr/bin/env python

"""
Run add gamma to restart file in order to run Hvrd's plant water stress function
"""

__author__    = "MU Mengyuan"

import os
import sys
import glob
import pandas as pd
import numpy as np
import netCDF4 as nc
import datetime
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.signal import savgol_filter

def main(file_type, restart_fname,para_fname,rate):

    # create file and write global attributes
    restart = nc.Dataset(restart_fname, "r+", format="NETCDF4")
    para    = nc.Dataset(para_fname, 'r')

    if file_type == "lis4real":

        # check relsmc
        # tmp = ( (restart.variables['SoilMoist'][3,:,:]-para.variables['swilt_inst'][0,:,:]) \
        #       / (para.variables['Porosity_inst'][0,:,:] - para.variables['swilt_inst'][0,:,:]) ).filled(-9999.) \
        #       - restart.variables['RelSMC'][3,:,:].filled(-9999.)
        # print(tmp.sum())


        # ignoring WATR (residual water content) since it is 0.0
        restart.variables['SoilMoist'][:,:,:] = \
               np.where( np.all([restart.variables['SoilMoist'][:,:,:]*rate < para.variables['Porosity_inst'][0:6,:,:],
                         restart.variables['SoilMoist'][:,:,:] != -9999.], axis=0),
                         restart.variables['SoilMoist'][:,:,:]*rate, para.variables['Porosity_inst'][0:6,:,:] )

        # ??? how to calculate it ???
        restart.variables['SoilWet'][:,:,:] = \
               np.where( np.all([restart.variables['SoilWet'][:,:,:]*rate < restart.variables['SoilMoist'][:,:,:],
                         restart.variables['SoilWet'][:,:,:] != -9999.], axis=0),
                         restart.variables['SoilWet'][:,:,:]*rate, restart.variables['SoilMoist'][:,:,:]    )


          avail_sm_mm = 0.
          do l=1,ms
             avail_smoist(l) = max(real(cable_struc(n)%cable(t)%wb(l) &
                  -cable_struc(n)%cable(t)%wbice(l)) - &
                  cable_struc(n)%cable(t)%swilt, 0.0)
             avail_sm_mm = avail_sm_mm + avail_smoist(l)*1000.0*&
                  cable_struc(n)%cable(t)%zse(l)
          enddo
          soilwet = avail_sm_mm/((cable_struc(n)%cable(t)%ssat-&
               cable_struc(n)%cable(t)%swilt)*&
               SUM(cable_struc(n)%cable(t)%zse(:)))



        restart.variables['RelSMC'][:,:,:] = ((restart.variables['SoilMoist'][:,:,:]-para.variables['swilt_inst'][0:6,:,:])\
                                              /(para.variables['Porosity_inst'][0:6,:,:] - para.variables['swilt_inst'][0:6,:,:])).filled(-9999.)


    elif file_type == "LIS_RST":

    restart.variables['WBTOT'][:] = np.where( np.all([restart.variables['SoilMoist'][:,:,:]*rate < para.variables['Porosity_inst'][0:6,:,:],
                     restart.variables['SoilMoist'][:,:,:] != -9999.], axis=0),
                     restart.variables['SoilMoist'][:,:,:]*rate, para.variables['Porosity_inst'][0:6,:,:] )

 WBTOT, WB, GWWB, WTD, WATR


    restart.close()
    para.close()

if __name__ == "__main__":

    file_type = "lis4real"
    rate      = 0.8

    if rate < 1.0:
        case_name = 'Princeton_ctl_watmove_new_veg_dry'
    elif rate >1.0:
        case_name = 'Princeton_ctl_watmove_new_veg_wet'

    if file_type == "lis4real":
        restart_fname = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/%s/bdy_data/lis4real_input.d01.nc" % case_name
    elif file_type == "LIS_RST"
        restart_fname = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/%s/LIS_output/LIS_RST_CABLE_200812310000.d01.nc" % case_name

    para_fname    = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/Princeton_ctl_para/LIS_output/LIS.CABLE.2008120100.d01.nc"


    main(file_type, restart_fname,para_fname, rate)
