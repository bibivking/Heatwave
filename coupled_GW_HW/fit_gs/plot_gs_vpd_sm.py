#!/usr/bin/env python

__author__  = "Mengyuan Mu"
__version__ = "2021-05-31"
__email__   = "mu.mengyuan815@gmail.com"

'''

'''

import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import datetime as dt
import netCDF4 as nc
from matplotlib import cm
from matplotlib import ticker
from scipy.interpolate import griddata
import scipy.stats as stats
from sklearn.metrics import mean_squared_error

def plot_gs_vpd_SM(fcable, ring):

    HPA_2_PA     = 100.0
    KPA_2_PA     = 1000.0
    DEG_2_KELVIN = 273.15
    SEC_2_D      = 24*60*60
    M_2_MM       = 1000.
    '''
    mol H20 m-2 s-1 -> mm s-1
    Tair = 20
    P    = 101.3
    '''
    molm2s1_2_mms = 18 * 0.000001

    # Read in data
    canht  = 20.
    press  = read_cable_var(fcable, "PSurf") # surface air pressure # hPa
    qair   = read_cable_var(fcable, "Qair")  # Surface specific humidity # kg/kg
    tair   = read_cable_var(fcable, "Tair")  # Surface air temperature # K
    rnet   = read_cable_var(fcable, "Rnet")  # Net radiation # W m-2
    wind   = read_cable_var(fcable, "Wind")  # Wind speed # m s-1
    # print(tair)

    trans  = read_obs_trans(ring)
    # print(trans)

    SM     = read_obs_swc_tdr(ring)
    # print(SM)

    press_tmp, trans_tmp= get_same_dates(press, trans)
    qair_tmp, trans_tmp = get_same_dates(qair, trans)
    tair_tmp, trans_tmp = get_same_dates(tair, trans)
    rnet_tmp, trans_tmp = get_same_dates(rnet, trans)
    wind_tmp, trans_tmp = get_same_dates(wind, trans)
    SM_tmp, trans_tmp   = get_same_dates(SM, trans)

    vpd    = qair_to_vpd(qair_tmp, tair_tmp, press_tmp)
    # print(vpd)
    # print(tair_tmp)

    press_tmp = press_tmp * HPA_2_PA
    tair_tmp  = tair_tmp - DEG_2_KELVIN
    vpd_tmp   = vpd * KPA_2_PA
    # print(trans_tmp)
    trans_tmp = trans_tmp / SEC_2_D / M_2_MM / molm2s1_2_mms
    # print(trans_tmp)

    # Coniferous forest, based on Jones 1976, from Jones '92, pg 67
    P  = PenmanMonteith(dz0v_dh=0.075, z0h_z0m=0.1)
    gc = np.zeros(np.size(vpd_tmp),dtype=np.float128)
    for i in np.arange(np.size(vpd_tmp)):
        # print(vpd_tmp[i], wind_tmp.values[i], rnet_tmp.values[i],
        #                  tair_tmp.values[i], press_tmp.values[i], trans_tmp.values[i])
        gc[i] = P.invert_penman(vpd_tmp[i], wind_tmp.values[i], rnet_tmp.values[i],
                         tair_tmp.values[i], press_tmp.values[i], trans_tmp.values[i], canht)
    print(gc)

    # SM <= 0.15
    gc1  =  gc[SM_tmp.values < 0.15]
    vpd1 =  vpd[SM_tmp.values < 0.15]

    # 0.15 <= SM_tmp.values < 0.3
    gc2  =  gc[np.all([SM_tmp.values < 0.3,SM_tmp.values >= 0.15],axis=0)]
    vpd2 =  vpd[np.all([SM_tmp.values < 0.3,SM_tmp.values >= 0.15],axis=0)]

    # 0.3 <= SM_tmp.values < 0.5
    gc3  =  gc[np.all([SM_tmp.values < 0.5,SM_tmp.values >= 0.3],axis=0)]
    vpd3 =  vpd[np.all([SM_tmp.values < 0.5,SM_tmp.values >= 0.3],axis=0)]

    # 0.5 <= SM_tmp.values < 0.7
    gc4  =  gc[np.all([SM_tmp.values < 0.7,SM_tmp.values >= 0.5],axis=0)]
    vpd4 =  vpd[np.all([SM_tmp.values < 0.7,SM_tmp.values >= 0.5],axis=0)]

    # 0.7 < SM_tmp.values <= 0.9
    gc5  =  gc[np.all([SM_tmp.values < 0.9,SM_tmp.values >= 0.7],axis=0)]
    vpd5 =  vpd[np.all([SM_tmp.values < 0.9,SM_tmp.values >= 0.7],axis=0)]

    # 0.9 < SM_tmp.values <= 1.
    gc6  =  gc[np.all([SM_tmp.values <= 1., SM_tmp.values >= 0.9],axis=0)]
    vpd6 =  vpd[np.all([SM_tmp.values <= 1., SM_tmp.values >= 0.9],axis=0)]

    # # fit solution
    # a, b = best_fit(Tobs['obs'][:], Tcable['cable'][:])
    # print(Tcable.values)

    x    = np.arange(0,7,0.01)
    # yfit = [a + b * xi for xi in x]

    # _____________ Make plot _____________
    fig = plt.figure(figsize=(6,5))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams["legend.markerscale"] = 2.0

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']  = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor']  = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    # set the box type of sequence number
    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
    # choose colormap
    colors = cm.Set2(np.arange(0,6))

    ax1 = fig.add_subplot(111)

    colors = cm.viridis(np.arange(0,6))

    # ax1.plot(x, x, lw=1.5, ls="-", c = almost_black, alpha=1.)
    # ax1.plot(x, yfit, lw=1.5, ls="--", c = almost_black, alpha=1.)
    ax1.scatter( vpd1, gc1, marker='o', c="darkslategray",edgecolors="darkslategray", s = 4., label="SM<0.15" , alpha=1)
    ax1.scatter( vpd2, gc2, marker='o', c="steelblue"    ,edgecolors="steelblue", s = 4., label="0.15<SM<0.3" , alpha=1)
    ax1.scatter( vpd3, gc3,marker='o',  c="limegreen"    ,edgecolors="limegreen", s = 4., label="0.3<SM<0.5" , alpha=1)
    ax1.scatter( vpd4, gc4, marker='o', c="springgreen"  ,edgecolors="springgreen", s = 4., label="0.5<SM<0.7", alpha=1)
    ax1.scatter( vpd5, gc5,marker='o',  c="gold"         ,edgecolors="gold", s = 4., label="0.7<SM<0.9" , alpha=1)
    # ax1.scatter( vpd6, gc6, marker='o', c="yellow"       ,edgecolors="yellow", s = 4., label="SM>0.9", alpha=0.5)


    ax1.axis('tight')
    # ax1.set_xticks([0,0.1,0.2,0.3,0.4])
    # ax1.set_xticklabels([0,0.1,0.2,0.3,0.4])
    # ax1.set_yticks([0,0.1,0.2,0.3,0.4])
    # ax1.set_yticklabels([0,0.1,0.2,0.3,0.4])
    ax1.set_xlim(0.,6.)
    # ax1.set_ylim(0.,0.38)
    ax1.set_xlabel('vpd (kPa)')
    ax1.set_ylabel('$g_{s}$ (-)')
    # ax1.text(0.82, 0.25, 'SM (-)', c = almost_black, transform=ax1.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    # 0.42
    ax1.legend(loc='upper right', frameon=False)

    fig.savefig("./plots/plot_gs_vpd_sm" , bbox_inches='tight', pad_inches=0.1)

def read_cable_var(fcable, var_name):

    """
    read a var from CABLE output
    """

    print("carry on read_cable_var")
    cable = nc.Dataset(fcable, 'r')
    Time  = nc.num2date(cable.variables['time'][:],cable.variables['time'].units)
    if var_name in ["TVeg", "ESoil", "Rainf", "GPP"]:
        var = pd.DataFrame(cable.variables[var_name][:,0,0]*1800., columns=['obs_a'])
    else:
        var = pd.DataFrame(cable.variables[var_name][:,0,0], columns=['obs_a'])
    var['Date'] = Time
    var = var.set_index('Date')
    if var_name in ["TVeg", "ESoil", "Rainf", "GPP"]:
        var = var.resample("D").agg('sum')
    else:
        print("is here")
        var = var.resample("D").agg('mean')
    var = var.sort_values(by=['Date'])
    var.index = np.arange(367,2924)

    return var

def read_obs_swc_tdr(ring):

    fobs   = "/srv/ccrc/data25/z5218916/data/Eucface_data/SM_2013-2019/eucSM1319_gap_filled.csv"
    tdr = pd.read_csv(fobs, usecols = ['Ring','Date','swc.tdr'])
    tdr['Date'] = pd.to_datetime(tdr['Date'],format="%d/%m/%Y",infer_datetime_format=False) # "%Y-%m-%d"
    tdr['Date'] = tdr['Date'] - pd.datetime(2011,12,31)
    tdr['Date'] = tdr['Date'].dt.days
    tdr = tdr.sort_values(by=['Date'])
    # divide neo into groups
    if ring == 'amb':
        subset = tdr[(tdr['Ring'].isin(['R2','R3','R6'])) & (tdr.Date > 366)]
    elif ring == 'ele':
        subset = tdr[(tdr['Ring'].isin(['R1','R4','R5'])) & (tdr.Date > 366)]
    else:
        subset = tdr[(tdr['Ring'].isin([ring]))  & (tdr.Date > 366)]

    subset = subset.groupby(by=["Date"]).mean()/100.
    subset = subset.rename({'swc.tdr' : 'obs_a'}, axis='columns')
    return subset

def read_obs_trans(ring):

    fobs_Trans = "/srv/ccrc/data25/z5218916/data/Eucface_data/FACE_PACKAGE_HYDROMET_GIMENO_20120430-20141115/data/Gimeno_wb_EucFACE_sapflow.csv"

    est_trans = pd.read_csv(fobs_Trans, usecols = ['Ring','Date','volRing'])
    est_trans['Date'] = pd.to_datetime(est_trans['Date'],format="%d/%m/%Y",infer_datetime_format=False)
    est_trans['Date'] = est_trans['Date'] - pd.datetime(2011,12,31)
    est_trans['Date'] = est_trans['Date'].dt.days
    est_trans = est_trans.sort_values(by=['Date'])
    # divide neo into groups
    if ring == 'amb':
       subs = est_trans[(est_trans['Ring'].isin(['R2','R3','R6'])) & (est_trans.Date > 366)]
    elif ring == 'ele':
       subs = est_trans[(est_trans['Ring'].isin(['R1','R4','R5'])) & (est_trans.Date > 366)]
    else:
       subs = est_trans[(est_trans['Ring'].isin([ring]))  & (est_trans.Date > 366)]

    subs = subs.groupby(by=["Date"]).mean()
    subs['volRing']   = subs['volRing'].clip(lower=0.)
    subs['volRing']   = subs['volRing'].replace(0., float('nan'))
    subs = subs.rename({'volRing' : 'obs_b'}, axis='columns')

    return subs

def get_same_dates(obs_a, obs_b):
    print("carry on get_same_dates")
    obs_a = obs_a['obs_a'].loc[obs_a.index.isin(obs_b.index)]
    obs_b = obs_b['obs_b'].loc[obs_b.index.isin(obs_a.index)]
    mask  = np.any([np.isnan(obs_a), np.isnan(obs_b)],axis=0)

    obs_a = obs_a[mask == False]
    obs_b = obs_b[mask == False]
    # print(obs_a, obs_b)

    return obs_a, obs_b

def best_fit(X, Y):

    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2
    # print(numer)
    # print(denum)
    b = numer / denum
    a = ybar - b * xbar
    # print(b)
    # print(a)
    print('best fit line: y = %2f + %2f x' %(a, b))

    return a, b

def qair_to_vpd(qair, tair, press):
    '''
    calculate vpd
    '''
    DEG_2_KELVIN = 273.15
    PA_TO_KPA = 0.001
    PA_TO_HPA = 0.01

    # convert back to Pa
    press_tmp = press / PA_TO_HPA
    tair_tmp  = tair - DEG_2_KELVIN

    # saturation vapor pressure
    es = 100.0 * 6.112 * np.exp((17.67 * tair_tmp) / (243.5 + tair_tmp))

    # vapor pressure
    ea = (qair * press_tmp) / (0.622 + (1.0 - 0.622) * tair_tmp)

    vpd = (es - ea) * PA_TO_KPA
    vpd = np.where(vpd < 0.05, 0.05, vpd)

    return vpd

class PenmanMonteith(object):

    """
    Penman-Monteith equation to calculate canopy transpiration.

    - Class also contain a method to invert canopy conducance (gc) if
      transpiration is already known

    References:
    -----------
    * Monteith and Unsworth (1990) Principles of Environmental
      Physics, pg. 247. Although I have removed the soil heat flux as G'DAY
      calculates soil evaporation seperately.

    """

    def __init__(self, dz0v_dh=0.075, z0h_z0m=0.1, d=0.67, use_ustar=False):

        """
        Parameters:
        -----------
        cp : float
            specific heat of dry air [MJ kg-1 degC-1]
        vk : float
            von Karman's constant [unitless]
        epsilon : float
            ratio molecular weight of water vap/dry air
        zele_sea : float
            elevation above sea level [m]
        dz0v_dh : float
            rate change of roughness for momentum with height
        displace_ratio : float
            zero plain displacement height
        z0h_z0m : float
            Ratio of the roughness length for heat to the roughness length for
            momentum, see comment in method below!!!
        """

        self.CP = 1010.0                 # specific heat of dry air (j kg-1 k-1)
        self.VK = 0.41                   # von Karman constan
        self.J_TO_MJ = 1.0E-6
        self.C_TO_K = 273.15
        self.dz0v_dh = dz0v_dh
        self.displace_ratio = d          # zero plan displacement height
        self.z0h_z0m = z0h_z0m
        self.RGAS = 8.314                # J mol-1 K-1
        self.H2OLV0 = 2.501E6            # latent heat H2O (J kg-1)
        self.H2OMW = 18E-3               # mol mass H20 (kg mol-1)
        self.MASS_AIR = 29.0E-3          # mol mass air (kg mol-1)
        self.use_ustar = use_ustar       # calc ga using ustar

    def calc_evaporation(self, vpd, wind, rnet, tair, press, gs, canht=None,
                         ustar=None, G=None):
        """
        Parameters:
        -----------
        vpd : float
            vapour pressure def [Pa]
        wind : float
            average daytime wind speed [m s-1]
        rnet : float
            net radiation [W m-2]
        tair : float
            temperature [degC]
        press : float
            average daytime pressure [Pa]
        gs : float
            stomatal conductance [mol m-2 s-1]
        canht : float
            canopy height [m]
        ustar : float
            friction velocity [m s-1]
        G : float
            soil heat flux [W m-2]

        Returns:
        --------
        trans : float
            transpiration [mol H20 m-2 s-1]
        """

        # use friction velocity
        if self.use_ustar:
            ga = self.calc_bdary_layer_conduc_from_ustar(wind, ustar, press,
                                                         tair)
        else:
            ga = self.canopy_bdary_layer_conduct(canht, wind, press, tair)
        lambdax = self.calc_latent_heat_of_vapourisation(tair)
        gamma = self.calc_pyschrometric_constant(lambdax, press)
        slope = self.calc_slope_of_sat_vapour_pressure_curve(tair)

        # Total leaf conductance to water vapour
        #gv = 1.0 / (1.0 / gs + 1.0 / ga)

        if G is None:
            G = 0.0 # ground heat flux

        arg1 = slope * (rnet - G) + ga * self.MASS_AIR * self.CP * vpd
        arg2 = slope + gamma * (1.0 + ga / gs)
        LE = arg1 / arg2 # W m-2

        # mol H20 m-2 s-1
        transpiration = LE / lambdax
        transpiration = np.where(transpiration < 0.0, 0.0, transpiration)

        return transpiration

    def invert_penman(self, vpd, wind, rnet, tair, press, trans, canht=None,
                      ustar=None, G=None):
        """
        Invert Penman-Monteith eqn to obtain canopy conductance

        Parameters:
        -----------
        vpd : float
            vapour pressure def [Pa]
        wind : float
            average daytime wind speed [m s-1]
        rnet : float
            net radiation [W m-2]
        tair : float
            temperature [degC]
        press : float
            average daytime pressure [Pa]
        trans : float
            transpiration [mol H20 m-2 s-1]
        canht : float
            canopy height [m]
        ustar : float
            friction velocity [m s-1]
        G : float
            soil heat flux [W m-2]

        Returns:
        --------
        gc : float
            canopy conductance [mol m-2 s-1]

        Reference:
        ---------
        * Landsberg and Sands, eqn 2.53
        """
        if self.use_ustar:
            ga = self.calc_bdary_layer_conduc_from_ustar(wind, ustar, press,
                                                         tair)
        else:
            ga = self.canopy_bdary_layer_conduct(canht, wind, press, tair)

        lambdax = self.calc_latent_heat_of_vapourisation(tair)
        gamma = self.calc_pyschrometric_constant(lambdax, press)
        slope = self.calc_slope_of_sat_vapour_pressure_curve(tair)
        lambda_E = trans * lambdax

        if G is None:
            G = 0.0 # ground heat flux

        arg1 = ga * gamma * lambda_E
        arg2 = slope * (rnet - G) - (slope + gamma) * lambda_E
        arg3 = ga * self.MASS_AIR * self.CP * vpd

        return arg1 / (arg2 + arg3)

    def calc_decoupling_coefficent(self, wind, tair, press, gs, canht=None,
                                  ustar=None):
        """
        Calculate decoupling coefficient.
            - As omega -> 0, leaf surface becomes strongly coupled to the
              atmosphere.
            - As omega -> 1, leaf surfaces are poorly coupled to the
              atmosphere.

        Parameters:
        -----------
        wind : float
            average daytime wind speed [m s-1]
        tair : float
            temperature [degC]
        press : float
            average daytime pressure [Pa]
        gs : float
            stomatal conductance [mol H20 m-2 s-1]
        canht : float
            canopy height [m]
        ustar : float
            friction velocity [m s-1]

        Returns:
        --------
        omega : float
            decoupling coefficient (-)


        References:
        -----------
        * McNaughton and Jarvis 1986
        """
        if self.use_ustar:
            ga = self.calc_bdary_layer_conduc_from_ustar(wind, ustar, press,
                                                         tair)
        else:
            ga = self.canopy_bdary_layer_conduct(canht, wind, press, tair)

        lambdax = self.calc_latent_heat_of_vapourisation(tair)
        gamma = self.calc_pyschrometric_constant(lambdax, press)
        slope = self.calc_slope_of_sat_vapour_pressure_curve(tair)

        epsilon = slope / gamma
        omega = (1.0 + epsilon) / (1.0 + epsilon + ga / gs)

        return (omega)

    def calc_bdary_layer_conduc_from_ustar(self, wind, ustar, press, tair):
        """
        Calculate boundary layer conductance using measured friction velocity,
        ustar

        Parameters:
        -----------
        wind : float
            average daytime wind speed [m s-1]
        ustar : float
            friction velocity [m s-1]

        Returns:
        --------
        ga : float
            canopy boundary layer conductance (mol m-2 s-1)

        References:
        ----------
        * Monteith & Unsworth p. 341 eqn. 17.8
        """

        # Convert from m s-1 to mol m-2 s-1
        # - note conversion in Jones '92 is mmol to mmol, but units cancel
        Tk = tair + self.C_TO_K
        cmolar = press / (self.RGAS * Tk)

        ga = 1.0 / (wind / ustar**2 + 6.2 * ustar**-0.667)
        ga *= cmolar

        return (ga)

    def canopy_bdary_layer_conduct(self, canht, wind, press, tair):
        """  Canopy boundary layer conductance, ga (from Jones 1992 p 68)

        Notes:
        ------
        'Estimates of ga for pine canopies from LAI of 3 to 6 vary from
        3.5 to 1.1 mol m-2 s-1  (Kelliher et al., 1993; Juang et al., 2007).'
        Drake et al, 2010, 17, pg. 1526.

        References:
        ------------
        * Jones 1992, pg. 67-8.
        * Monteith and Unsworth (1990), pg. 248. Note this in the inverted form
          of what is in Monteith (ga = 1 / ra)
        * Allen et al. (1989) pg. 651.
        * Gash et al. (1999) Ag forest met, 94, 149-158.

        Parameters:
        -----------
        wind : float
            average daytime wind speed [m s-1]
        press : float
            atmospheric pressure (Pa)
        tair : float
            air temperature (deg C)
        canht : float
            canopy height (m)

        Returns:
        --------
        ga : float
            canopy boundary layer conductance [mol m-2 s-1]
        """

        # Convert from m s-1 to mol m-2 s-1
        # - note conversion in Jones '92 is mmol to mmol, but units cancel
        Tk = tair + self.C_TO_K
        cmolar = press / (self.RGAS * Tk)

        # roughness length for momentum
        z0m = self.dz0v_dh * canht

        # z0h roughness length governing transfer of heat and vapour [m]
        # *Heat tranfer typically less efficent than momentum transfer. There is
        #  a lot of variability in values quoted for the ratio of these two...
        #  JULES uses 0.1, Campbell and Norman '98 say z0h = z0m / 5. Garratt
        #  and Hicks, 1973/ Stewart et al '94 say z0h = z0m / 7. Therefore for
        #  the default I am following Monteith and Unsworth, by setting the
        #  ratio to be 1, the code below is identical to that on page 249,
        #  eqn 15.7
        z0h = self.z0h_z0m * z0m

        # zero plan displacement height [m]
        d = self.displace_ratio * canht

        arg1 = self.VK**2 * wind
        arg2 = np.log((canht - d) / z0m)
        arg3 = np.log((canht - d) / z0h)

        ga = (arg1 / (arg2 * arg3)) * cmolar

        return (ga)

    def calc_slope_of_sat_vapour_pressure_curve(self, tair):
        """ Constant slope in Penman-Monteith equation

        Parameters:
        -----------
        tavg : float
            average daytime temperature

        Returns:
        --------
        slope : float
            slope of saturation vapour pressure curve [Pa K-1]
        """
        arg1 = self.calc_sat_water_vapour_press(tair+0.1)
        arg2 = self.calc_sat_water_vapour_press(tair)
        slope = (arg1 - arg2) / 0.1

        return (slope)

    def calc_sat_water_vapour_press(self, tac):
        """ Calculate saturated water vapour pressure at temperature TAC

        Parameters:
        -----------
        tac : float
            Celsius

        Returns:
        --------
        sat : float
            units: Pa

        References:
        -----------
        * Jones 1992 p 110 (note error in a - wrong units)
        """

        return (613.75 * np.exp(17.502 * tac / (240.97 + tac)));

    def calc_pyschrometric_constant(self, lambdax, press):
        """ Psychrometric constant ratio of specific heat of moist air at
        a constant pressure to latent heat of vaporisation.

        Parameters:
        -----------
        press : float
            air pressure (Pa)
        lambda : float
             latent heat of water vaporization (J mol-1)

        Returns:
        --------
        gamma : float
            pyschrometric constant [Pa K-1]
        """

        return (self.CP * self.MASS_AIR * press / lambdax)

    def calc_latent_heat_of_vapourisation(self, tair):
        """
        Latent heat of water vapour at air temperature

        Returns:
        -----------
        lambda : float
            latent heat of water vaporization [J mol-1]
        """

        return ((self.H2OLV0 - 2.365E3 * tair) * self.H2OMW)


if __name__ == "__main__":

    ring    = "amb"#"amb"
    fcable  = "/srv/ccrc/data25/z5218916/cable/EucFACE/EucFACE_run/outputs/met_LAI-08_vrt_swilt-watr-ssat_hyds10_31uni_teuc_sres_watr/EucFACE_amb_out.nc"

    plot_gs_vpd_SM(fcable, ring)
