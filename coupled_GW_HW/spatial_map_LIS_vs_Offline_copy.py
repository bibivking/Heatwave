#!/usr/bin/python

from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_scale_offline
from common_utils import get_reverse_colormap

def tree_mask(file_path,pft_var_name):

    var = Dataset(file_path, mode='r')
    if pft_var_name == "Landcover_inst":
        pft = var.variables[pft_var_name][0,:,:]
    elif pft_var_name == "iveg":
        pft = var.variables[pft_var_name][:,:]

    return pft

def mask_by_lat_lon(file_path, loc_lat, loc_lon, lat_name, lon_name):

    file = nc.Dataset(file_path, mode='r')
    lat  = file.variables[lat_name][:]
    lon  = file.variables[lon_name][:]
    print(lat)
    print(lon)
    if len(np.shape(lat)) == 1:
        print("len(np.shape(lat)) == 1")
        lat_spc = lat[1] - lat[0]
        lon_spc = lon[1] - lon[0]
        lons, lats = np.meshgrid(lon, lat)
        mask  = (lats > (loc_lat - lat_spc/2)) & (lats < (loc_lat + lat_spc/2)) & (lons > (loc_lon - lon_spc/2)) & (lons < (loc_lon + lon_spc/2))
    elif len(np.shape(lat)) == 2:
        print("len(np.shape(lat)) == 2")
        ### caution: lat=100, lon=100 is a random pixel, lis run over a small domain may not have such a point
        lat_spc = lat[100,100] - lat[99,100]
        lon_spc = lon[100,100] - lon[100,99]
        print(lat_spc)
        print(lon_spc)
        ### caution: due to irregular space in lis, using lat/lon +lat/lon_spc/2 may includes more than 1 pixel.
        ### I therefore let the space divied by 2.1 rather than 2
        mask  = (lat > (loc_lat - lat_spc/2.1)) & (lat < (loc_lat + lat_spc/2.1)) & (lon > (loc_lon - lon_spc/2.1)) & (lon < (loc_lon + lon_spc/2.1))
    return mask

def plot_map_var_offline(file_path, var_names, is_lnd=False, layer=None, is_tree=False):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var = Dataset(file_path, mode='r')

    time = var.variables['time']
    lons = var.variables['longitude'][:,:] # or ['lon']
    lats = var.variables['latitude'][:,:]  # or ['lat']
    # lon, lat = np.meshgrid(lons, lats)
    print(lons)
    print(lats)

    for var_name in var_names:
        print(var_name)
        scale, units = get_land_var_scale_offline(var_name)

        Var     = np.mean(var.variables[var_name],axis=0)*scale
        def_val = var.variables[var_name]._FillValue
        if is_tree:
            pft = tree_mask(file_path,"iveg")
            Var = np.where(pft <=4, Var, def_val)

        print(Var.shape)

        fig = plt.figure(figsize=(7,5))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([140,154,-40,-28])

        ax.coastlines(resolution="50m",linewidth=1)
        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator([140,145,150])
        gl.ylocator = mticker.FixedLocator([-40,-35,-30])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black'}
        gl.ylabel_style = {'size':10, 'color':'black'}
        # Plot windspeed

        if is_lnd:
            print("is_lnd"+str(is_lnd))
            if var_name in ["iveg","isoil"]:
                clevs = np.linspace(0.5,15.5, num=16)
            elif var_name in ["Albedo","sfc","ssat","swilt"]:
                clevs = np.linspace(0.0,0.5, num=21)
            elif var_name in ["sand","clay","silt"]:
                clevs = np.linspace(0.0,1.0, num=21)
            elif var_name in ["LAI"]:
                clevs = np.linspace(0.,8., num=16)
            elif var_name in ["bch"]:
                clevs = np.linspace(0.,13., num=26)
            elif var_name in ["hyds"]:
                clevs = np.linspace(0.,0.05, num=21)
            elif var_name in ["sucs"]:
                clevs = np.linspace(0.,1000., num=21)
            else:
                clevs = np.linspace(np.min(Var),np.max(Var), num=21)
        else:
            if var_name == 'WatTable':
                clevs = np.arange(0,15,1)
            elif var_name == 'SoilMoist':
                clevs = np.arange(0,0.5,0.05)
            else:
                clevs = np.arange(0,5.,0.5)

        if len(np.shape(Var)) == 2:
            plt.contourf(lons, lats, Var, clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG)#, extend='both')#seismic) # T2M_daily_avg
        else:
            plt.contourf(lons, lats, Var[layer,:,:], clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG)#, extend='both')#seismic) # T2M_daily_avg

        plt.title(var_name, size=16)
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.set_label(units,size=14,rotation=270,labelpad=15)
        cb.ax.tick_params(labelsize=10)

        if is_tree:
            message = var_name+"_tree"
        else:
            message = var_name

        if len(np.shape(Var)) == 2:
            plt.savefig('./plots/lis_vs_offline_26Sep/spatial_map_offline_'+message+'.png',dpi=300)
        else:
            plt.savefig('./plots/lis_vs_offline_26Sep/spatial_map_offline_'+message+'_lvl-'+str(layer)+'.png',dpi=300)

def plot_map_var_lis(file_path, wrf_path, case_name, var_names, is_lnd=False, layer=None, is_tree=False):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:

    var = Dataset(file_path, mode='r')

    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    for var_name in var_names:

        print(var_name)
        scale, units = get_land_var_scale(var_name)

        Var   = np.mean(var.variables[var_name],axis=0)*scale
        print(Var.shape)

        if is_tree:
            pft = tree_mask(file_path,"Landcover_inst")
            Var = np.where(pft <= 4, Var, -9999.)

        # Make plots
        fig = plt.figure(figsize=(7,5))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([140,154,-40,-28])
        ax.coastlines(resolution="50m",linewidth=1)

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator([140,145,150])
        gl.ylocator = mticker.FixedLocator([-40,-35,-30])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black'}
        gl.ylabel_style = {'size':10, 'color':'black'}

        if is_lnd:
            print("is_lnd"+str(is_lnd))
            if var_name in ["Landcover_inst","Soiltype_inst"]:
                clevs = np.linspace(0.5,15.5, num=16)
            elif var_name in ["Albedo_inst","SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst"]:
                clevs = np.linspace(0.0,0.5, num=21)
            elif var_name in ["SandFrac_inst","ClayFrac_inst","SiltFrac_inst"]:
                clevs = np.linspace(0.0,1.0, num=21)
            elif var_name in ["LAI_inst"]:
                clevs = np.linspace(0.,7.5, num=17)
            elif var_name in ["Bch_inst"]:
                clevs = np.linspace(0.,13., num=27)
            elif var_name in ["Hyds_inst"]:
                clevs = np.linspace(0.,0.05, num=21)
            elif var_name in ["Sucs_inst"]:
                clevs = np.linspace(0.,1000., num=21)
            else:
                clevs = np.linspace(np.min(Var),np.max(Var), num=21)
        else:
            if var_name == 'WaterTableD_tavg':
                clevs = np.arange(0,15,1)
            elif var_name == 'SoilMoist_inst':
                clevs = np.arange(0,0.5,0.05)
            else:
                clevs = np.arange(0,5,0.5)

        print(Var.shape)
        if len(np.shape(Var)) == 2:
            plt.contourf(lon, lat, Var[:,:], clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG)#seismic)
        else:
            plt.contourf(lon, lat, Var[layer,:,:], clevs, transform=ccrs.PlateCarree(),cmap=plt.cm.BrBG)#seismic)
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

        if units == None:
            units_string = var.variables[var_name].units
        else:
            units_string = units

        plt.title(var_name+' ('+units_string+')', size=16)

        cb.ax.tick_params(labelsize=10)
        cb.set_label(units_string,size=14,rotation=270,labelpad=15)

        if is_tree:
            message = var_name+'_'+case_name+"_tree"
        else:
            message = var_name+'_'+case_name

        if len(np.shape(Var)) == 2:
            plt.savefig('./plots/lis_vs_offline_26Sep/spatial_map_lis_'+message+'.png',dpi=300)
        else:
            plt.savefig('./plots/lis_vs_offline_26Sep/spatial_map_lis_'+message+'_lvl-'+str(layer)+'.png',dpi=300)

        Var = None

def statistic_bar(file_path_offline, file_path_lis, wrf_path, is_tree=False, loc_lat=None, loc_lon=None):

    off_stats   = np.zeros(7)
    lis_stats   = np.zeros(7)
    scale       = 3600*24.*366.

    # ================= Read offline ==================
    offline = Dataset(file_path_offline, mode='r')

    rain_off = offline.variables['Rainf'][0,:,:]
    if loc_lat ==None:
        if is_tree:
            pft = tree_mask(file_path_offline,"iveg")
            mask_off = (rain_off >= 0.) & (pft <= 4)
        else:
            mask_off = (rain_off >= 0.)

        rain_off_tmp = np.mean(offline.variables['Rainf'],axis=0)

        off_stats[0] = np.nanmean(np.where(rain_off_tmp[mask_off], rain_off_tmp[mask_off], np.nan))

        evap_off_tmp = np.mean(offline.variables['Evap'],axis=0)
        off_stats[1] = np.nanmean(np.where(evap_off_tmp[mask_off], evap_off_tmp[mask_off], np.nan))

        tveg_off_tmp = np.mean(offline.variables['TVeg'],axis=0)
        off_stats[2] = np.nanmean(np.where(tveg_off_tmp[mask_off], tveg_off_tmp[mask_off], np.nan))

        esoil_off_tmp = np.mean(offline.variables['ESoil'],axis=0)
        off_stats[3]  = np.nanmean(np.where(esoil_off_tmp[mask_off], esoil_off_tmp[mask_off], np.nan))

        ecanop_off_tmp = np.mean(offline.variables['ECanop'],axis=0)
        off_stats[4]   = np.nanmean(np.where(ecanop_off_tmp[mask_off], ecanop_off_tmp[mask_off], np.nan))

        qs_off_tmp     = np.mean(offline.variables['Qs'],axis=0)
        off_stats[5]   = np.nanmean(np.where(qs_off_tmp[mask_off], qs_off_tmp[mask_off], np.nan))

        qsb_off_tmp    = np.mean(offline.variables['Qsb'],axis=0)
        off_stats[6]   = np.nanmean(np.where(qsb_off_tmp[mask_off], qsb_off_tmp[mask_off], np.nan))
    else:
        print(offline.variables['Rainf'])
        off_stats[0] = np.mean(offline.variables['Rainf'],axis=0)
        off_stats[1] = np.mean(offline.variables['Evap'],axis=0)
        off_stats[2] = np.mean(offline.variables['TVeg'],axis=0)
        off_stats[3] = np.mean(offline.variables['ESoil'],axis=0)
        off_stats[4] = np.mean(offline.variables['ECanop'],axis=0)
        off_stats[5] = np.mean(offline.variables['Qs'],axis=0)
        off_stats[6] = np.mean(offline.variables['Qsb'],axis=0)

    # ================= Read lis-cable =================
    lis     = Dataset(file_path_lis, mode='r')
    wrf     = Dataset(wrf_path,  mode='r')

    # mask: keep SE Aus
    lon     = wrf.variables['XLONG'][0,:,:]
    lat     = wrf.variables['XLAT'][0,:,:]
    rain    = lis.variables['Rainf_tavg'][0,:,:]
    if loc_lat ==None:
        if is_tree:
            pft = tree_mask(file_path_lis,"Landcover_inst")
            mask    = (lon >= 140) & (lon <= 154) & (lat >= -40) & (lat <= -28) & (rain >= 0.) & (pft <= 4)
        else:
            mask    = (lon >= 140) & (lon <= 154) & (lat >= -40) & (lat <= -28) & (rain >= 0.)

        rain_lis_tmp = np.mean(lis.variables['Rainf_tavg'],axis=0)
        lis_stats[0] = np.nanmean(np.where(rain_lis_tmp[mask], rain_lis_tmp[mask], np.nan))

        evap_lis_tmp = np.mean(lis.variables['Evap_tavg'],axis=0)
        lis_stats[1] = np.nanmean(np.where(evap_lis_tmp[mask], evap_lis_tmp[mask], np.nan))

        tveg_lis_tmp = np.mean(lis.variables['TVeg_tavg'],axis=0)
        lis_stats[2] = np.nanmean(np.where(tveg_lis_tmp[mask], tveg_lis_tmp[mask], np.nan))

        esoil_lis_tmp = np.mean(lis.variables['ESoil_tavg'],axis=0)
        lis_stats[3]  = np.nanmean(np.where(esoil_lis_tmp[mask], esoil_lis_tmp[mask], np.nan))

        ecanop_lis_tmp = np.mean(lis.variables['ECanop_tavg'],axis=0)
        lis_stats[4]   = np.nanmean(np.where(ecanop_lis_tmp[mask], ecanop_lis_tmp[mask], np.nan))

        qs_lis_tmp     = np.mean(lis.variables['Qs_tavg'],axis=0)
        lis_stats[5]   = np.nanmean(np.where(qs_lis_tmp[mask], qs_lis_tmp[mask], np.nan))

        qsb_lis_tmp    = np.mean(lis.variables['Qsb_tavg'],axis=0)
        lis_stats[6]   = np.nanmean(np.where(qsb_lis_tmp[mask], qsb_lis_tmp[mask], np.nan))

    else:
        mask     = mask_by_lat_lon(file_path_lis, loc_lat, loc_lon, "lat", "lon")
        print(mask)
        rain_lis_tmp = np.mean(lis.variables['Rainf_tavg'],axis=0)
        print(rain_lis_tmp)
        lis_stats[0] = rain_lis_tmp[mask]

        evap_lis_tmp = np.mean(lis.variables['Evap_tavg'],axis=0)
        lis_stats[1] = evap_lis_tmp[mask]

        tveg_lis_tmp = np.mean(lis.variables['TVeg_tavg'],axis=0)
        lis_stats[2] = tveg_lis_tmp[mask]

        esoil_lis_tmp = np.mean(lis.variables['ESoil_tavg'],axis=0)
        lis_stats[3]  = esoil_lis_tmp[mask]

        ecanop_lis_tmp = np.mean(lis.variables['ECanop_tavg'],axis=0)
        lis_stats[4]   = ecanop_lis_tmp[mask]

        qs_lis_tmp     = np.mean(lis.variables['Qs_tavg'],axis=0)
        lis_stats[5]   = qs_lis_tmp[mask]

        qsb_lis_tmp    = np.mean(lis.variables['Qsb_tavg'],axis=0)
        lis_stats[6]   = qsb_lis_tmp[mask]

    off_stats = off_stats*scale
    lis_stats = lis_stats*scale

    # plotting
    labels = ['Rain', 'Evap', 'TVeg', 'Esoil', 'Ecan','Qs','Qsb']

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, off_stats, width, label='Off')
    rects2 = ax.bar(x + width/2, lis_stats, width, label='LIS')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('mm/yr')
    ax.set_title('Water Balance')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    # ax.bar_label(rects1, padding=3)
    # ax.bar_label(rects2, padding=3)

    fig.tight_layout()
    message = ""

    if is_tree:
        message = message + "_tree"
    if loc_lat != None:
        message = message + "_lat="+str(loc_lat) + "_lon="+str(loc_lon)

    plt.savefig('./plots/lis_vs_offline_26Sep/water_balance_lis_vs_off_statistic'+message+'.png',dpi=300)

def time_series_single_pixel(file_path_offline, file_path_lis, varname_off, varname_lis, loc_lat=None, loc_lon=None):

    offline    = Dataset(file_path_offline, mode='r')
    times_off  = offline.variables['time']
    ntime_off  = len(times_off)

    lis        = Dataset(file_path_lis, mode='r')
    times_lis  = offline.variables['time']
    ntime_lis  = len(times_lis)
    var_lis    = np.zeros(ntime_lis)
    # print(np.shape(var_lis))
    mask     = mask_by_lat_lon(file_path_lis, loc_lat, loc_lon, "lat", "lon")

    for i in np.arange(len(varname_off)):
        print(varname_off[i])
        scale, units = get_land_var_scale_offline(varname_off[i])
        if varname_off[i] == "SoilMoist":
            for l in np.arange(6):
                var_off   = offline.variables[varname_off[i]][:,l,0,0]*scale
                for j in np.arange(ntime_lis):
                    scale, units = get_land_var_scale(varname_lis[i])
                    var_lis[j] = np.mean(lis.variables[varname_lis[i]][j,l][mask])*scale

                x1 = np.arange(ntime_off)
                x2 = np.arange(ntime_lis)

                fig, ax = plt.subplots()
                ax.plot(x1, var_off, c = "blue",   label='Off',alpha=0.5)
                ax.plot(x2, var_lis, c = "orange", label='Lis',alpha=0.5)

                # ax.set_ylabel('mm')
                ax.set_title(varname_off[i])
                ax.set_xticks(x1[::1440])
                ax.set_xticklabels(np.arange(1,366,30))
                ax.legend()

                fig.tight_layout()
                message = varname_off[i]+"lvl="+str(l)
                if loc_lat != None:
                    message = message + "_lat="+str(loc_lat) + "_lon="+str(loc_lon)

                plt.savefig('./plots/lis_vs_offline_26Sep/time_series_lis_vs_off_'+message+'.png',dpi=300)

        else:
            var_off   = offline.variables[varname_off[i]][:,0,0]*scale
            for j in np.arange(ntime_lis):
                scale, units = get_land_var_scale(varname_lis[i])
                var_lis[j] = np.mean(lis.variables[varname_lis[i]][j][mask])*scale

            x1 = np.arange(ntime_off)
            x2 = np.arange(ntime_lis)

            fig, ax = plt.subplots()
            ax.plot(x1, var_off, c = "blue",   label='Off',alpha=0.5)
            ax.plot(x2, var_lis, c = "orange", label='Lis',alpha=0.5)

            # ax.set_ylabel('mm')
            ax.set_title(varname_off[i])
            ax.set_xticks(x1[::1440])
            ax.set_xticklabels(np.arange(1,366,30))
            ax.legend()

            fig.tight_layout()
            message = varname_off[i]
            if loc_lat != None:
                message = message + "_lat="+str(loc_lat) + "_lon="+str(loc_lon)

            plt.savefig('./plots/lis_vs_offline_26Sep/time_series_lis_vs_off_'+message+'.png',dpi=300)

def print_compare_land_info(file_path_offline, file_path_lis, varname_off, varname_lis, loc_lat=None, loc_lon=None):

    offline    = Dataset(file_path_offline, mode='r')
    lis        = Dataset(file_path_lis, mode='r')
    mask     = mask_by_lat_lon(file_path_lis, loc_lat, loc_lon, "lat", "lon")

    for i in np.arange(len(varname_off)):
        print(varname_off[i])
        if varname_off[i] in ["sand","clay","silt","sfc","ssat","swilt","hyds","bch","sucs"]:
            var_off   = offline.variables[varname_off[i]][:,0,0]
            var_lis   = lis.variables[varname_lis[i]][0][mask]

            print("off is")
            print(var_off)
            print("lis is")
            print(var_lis)
        elif varname_off[i] in ["iveg","isoil",'elev',"Albedo"]:
            var_off   = offline.variables[varname_off[i]][0,0]
            var_lis   = lis.variables[varname_lis[i]][0][mask]

            print("off is")
            print(var_off)
            print("lis is")
            print(var_lis)
        elif varname_off[i] == "LAI":
            var_lis   = np.zeros(17568)
            var_off   = offline.variables[varname_off[i]][:,0,0]
            for j in np.arange(17568):
                var_lis[j]   = lis.variables[varname_lis[i]][j][mask]

            print("off-lis is")
            print(np.sum(var_off - var_lis))

        var_off = None
        var_lis = None

if __name__ == "__main__":

    var_LIS_names      = ['FWsoil_tavg','Qs_tavg','Qsb_tavg','WaterTableD_tavg','Qle_tavg','Qh_tavg','Qg_tavg',
                          'GWwb_tavg','Rainf_tavg','Evap_tavg','ESoil_tavg','ECanop_tavg','TVeg_tavg']

    var_offline_names  = ['Fwsoil','Qs','Qsb','WatTable','Qle','Qh','Qg','GWMoist','Rainf','Evap','ESoil','ECanop','TVeg']

    var_LIS_met_names  = ['Psurf_f_inst'] #'Rainf_f_inst','Snowf_tavg','Qair_f_inst','Wind_f_inst','SWdown_f_inst','LWdown_f_inst','Tair_f_inst'

    var_offline_met_names  = ['PSurf'] #'Rainf','Snowf','Qair','Wind','SWdown','LWdown','Tair'

    var_LIS_soil_names     = ["SoilMoist_inst"]

    var_offline_soil_names = ["SoilMoist"]

    var_LIS_landinfo_names =  [ "SandFrac_inst","ClayFrac_inst","SiltFrac_inst",
                                "SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst","Hyds_inst","Bch_inst","Sucs_inst",
                                "LAI_inst","Landcover_inst","Soiltype_inst","Elevation_inst" ,"Albedo_inst"]

    var_offline_landinfo_names = ["sand","clay","silt",
                                  "sfc","ssat","swilt","hyds","bch","sucs",
                                  "LAI","iveg","isoil",'elev',"Albedo"]

    ##########################
    #   Plot offline CABLE   #
    # ##########################
    is_tree   = False
    is_lnd    = False
    path      = '/g/data/w35/mm3972/model/cable/runs/test_para_chg_dpt/uniform_6layer_fix_satfrac_mmy_trunk/outputs/'
    file_name = 'cable_out_2000_SE_Aus.nc'
    file_path = path + file_name
    var_names = var_offline_names #var_offline_soil_names #var_offline_names #var_offline_landinfo_names #var_offline_soil_names # var_offline_names # var_offline_soil_names
    plot_map_var_offline(file_path, var_names, is_lnd=is_lnd, is_tree=is_tree )
    var_names = var_offline_soil_names
    for layer in np.arange(6):
       plot_map_var_offline(file_path, var_names, is_lnd=is_lnd, is_tree=is_tree, layer=layer)
    layer = None
    plot_map_var_offline(file_path, var_offline_names, is_lnd=is_lnd, is_tree=is_tree, layer=layer)

    ##############################
    #   plot plot_map_var_lis    #
    ##############################
    is_tree    = False
    is_lnd     = False
    case_name  = 'ctl_26Sep'
    file_name  = "LIS.CABLE.200001-200012.nc"
    layer      = None
    path       = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_26Sep/LIS_output/"
    file_path  = path + file_name
    wrf_path   = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_14Aug/WRF_output_copy/wrfout_d01_2013-01-01_03:00:00"
    # var_names  = var_LIS_soil_names #var_LIS_names #var_LIS_landinfo_names #var_LIS_soil_names # var_LIS_names # var_LIS_soil_names
    # # plot_map_var_lis(file_path, wrf_path, case_name, var_names, is_lnd=is_lnd, is_tree=is_tree)
    # for layer in np.arange(6):
    #     plot_map_var_lis(file_path, wrf_path, case_name, var_names, is_lnd=is_lnd, is_tree=is_tree, layer=layer)
    plot_map_var_lis(file_path, wrf_path, case_name, var_LIS_names, is_lnd=is_lnd, is_tree=is_tree, layer=layer)

    ##############################
    #   plot statistic_bar    #
    ##############################
    is_tree           = False

    loc_lat           = -34 #-24.255707
    loc_lon           = 145  #135.95001

    file_path_offline = '/g/data/w35/mm3972/model/cable/runs/test_para_chg_dpt/uniform_6layer_fix_satfrac_mmy_trunk/outputs/cable_out_2000_SE_Aus.nc'
    # file_path_offline = '/g/data/w35/mm3972/model/cable/runs/pixel_comp_lis/outputs/ERAI_05hr_pixel_met_from_LIS-CABLE_satfrac_fixed_output.nc'
    file_path_lis     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_26Sep/LIS_output/LIS.CABLE.200001-200012.nc'
    wrf_path          = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_14Aug/WRF_output_copy/wrfout_d01_2013-01-01_03:00:00"
    statistic_bar(file_path_offline, file_path_lis, wrf_path, is_tree=is_tree)#, loc_lat=loc_lat, loc_lon=loc_lon)


    # ########################################
    # #   plot time_series_single_pixel      #
    # ########################################
    # file_path_offline = '/g/data/w35/mm3972/model/cable/runs/pixel_comp_lis/outputs/ERAI_05hr_pixel_met_from_LIS-CABLE_satfrac_fixed_output.nc'
    # file_path_lis     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_26Sep/LIS_output/LIS.CABLE.200001-200012.nc'
    # loc_lat           = -34 #-24.255707
    # loc_lon           = 145  #135.95001
    # time_series_single_pixel(file_path_offline, file_path_lis, var_offline_names, var_LIS_names, loc_lat=loc_lat, loc_lon=loc_lon)


    # # ########################################
    # # #   compare land parameters            #
    # # ########################################
    # file_path_offline = '/g/data/w35/mm3972/model/cable/runs/pixel_comp_lis/outputs/ERAI_05hr_pixel_met_from_LIS-CABLE_satfrac_fixed_output.nc'
    # file_path_lis     = '/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ctl_15Sep/LIS_output_large_domain_half-hourly/LIS.CABLE.200001-200012.nc'
    # loc_lat           = -34 #-24.255707
    # loc_lon           = 145  #135.95001
    # print_compare_land_info(file_path_offline, file_path_lis, var_offline_landinfo_names, var_LIS_landinfo_names, loc_lat=loc_lat, loc_lon=loc_lon)
