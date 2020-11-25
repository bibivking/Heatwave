;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

undef("regrid_data")
function regrid_data(var,interp_method,lat2d,lon2d,lat,lon,src_mask,dst_mask)
local var_regrid
begin

    Opt                   = True

    Opt@SrcFileName       = "WRF_SCRIP.nc"      ; Name of source and ???
    Opt@DstFileName       = "1degree.nc"        ; destination files  ???

    Opt@InterpMethod      = interp_method

    Opt@SrcRegional       = True
    Opt@SrcGridMask       = src_mask

    Opt@DstRegional       = True
    ; Opt@DstLLCorner       = (/ minlat,minlon /)
    ; Opt@DstURCorner       = (/ maxlat,maxlon /)
    Opt@DstGridType       = "rectilinear"
    Opt@DstGridLat        = lat
    Opt@DstGridLon        = lon
    Opt@DstGridMask       = dst_mask

    Opt@ForceOverwrite    = True
    Opt@PrintTimings      = True
    Opt@Debug             = True

    var@lat2d = lat2d              ; This information will be used by
    var@lon2d = lon2d              ; ESMF_regrid for the source grid
    var_regrid = ESMF_regrid(var,Opt)     ; Do the regridding for TMP
    printVarSummary(var_regrid)
    return(var_regrid)
end

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

; =================================== Path =====================================
case_name = "CTL_vs_PUMP"
path1 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_GSWP3_GW/outputs/CTL"
path2 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_GSWP3_GW_pumping/outputs/PUMP"

path_GLEAM= "/g/data/w35/Shared_data/Observations/Global_ET_products/GLEAM_v3_3/3_3a/yearly";
path_CLASS= "/g/data/w35/mm3972/data/CLASS_v1"
path_DOLCE= "/g/data/w35/mm3972/data/DOLCE" ; DOLCE_v1.0_2000.nc
path_LORA = "/g/data/w35/mm3972/data/LORA"  ; LORA_v1.0_2000.nc
path_CERES= "/g/data/w35/mm3972/data/CERES" ; CERES_EBAF-Surface_Ed2.8_Subset_200003-201502.nc

; ================================= Predefine ==================================
year_s = 2003
year_e = 2008
year_sum = year_e-year_s+1

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Soil_thickness1 = Soil_thickness
;(/0.005, 0.075, 0.154, 0.409, 1.085, 2.872/) ; LIS-CABLE


; ================================== OPTIONS ===================================
plot_type = "all"
scale = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"

filename = path1+"/cable_out_"+year_s+".nc"
f = addfile (filename,"r")

lat_size  = dimsizes(f->y)
lon_size  = dimsizes(f->x)
CABLE_lat = f->y
CABLE_lon = f->x
default_value = -1.e+33

if scale .eq. "Global" then
  A = -90.
  B = 90.
  C = 0.
  D = 360.
else if scale .eq. "SE-AU" then
  A = -40.
  B = -27.
  C = 140.
  D = 155.
else if scale .eq. "AU" then
  A = -44.
  B = -10.
  C = 112.
  D = 154.
end if
end if
end if
delete(f)


filename1 = path1+"/cable_out_"+year_s+".nc"
f1        = addfile (filename1,"r")

Rain_CTL     = dim_sum_n_Wrap(f1->Rainf(:,{A:B},{C:D}),0)
Tair_CTL     = dim_avg_n_Wrap(f1->Tair(:,{A:B},{C:D}),0)
Qair_CTL     = dim_avg_n_Wrap(f1->Qair(:,{A:B},{C:D}),0)
Psurf_CTL    = dim_avg_n_Wrap(f1->PSurf(:,{A:B},{C:D}),0)
SWdown_CTL   = dim_avg_n_Wrap(f1->SWdown(:,{A:B},{C:D}),0)
LWdown_CTL   = dim_avg_n_Wrap(f1->LWdown(:,{A:B},{C:D}),0)
Wind_CTL     = dim_avg_n_Wrap(f1->Wind(:,{A:B},{C:D}),0)

Evap_CTL     = dim_sum_n_Wrap(f1->Evap(:,{A:B},{C:D}),0)
TVeg_CTL     = dim_sum_n_Wrap(f1->TVeg(:,{A:B},{C:D}),0)
ESoil_CTL    = dim_sum_n_Wrap(f1->ESoil(:,{A:B},{C:D}),0)
ECanop_CTL   = dim_sum_n_Wrap(f1->ECanop(:,{A:B},{C:D}),0)
Runoff_CTL   = dim_sum_n_Wrap(f1->Qs(:,{A:B},{C:D}),0)
Qrec_CTL     = dim_sum_n_Wrap(f1->Qrecharge(:,{A:B},{C:D}),0)
Rnet_CTL     = dim_avg_n_Wrap(f1->Rnet(:,{A:B},{C:D}),0)
Qle_CTL      = dim_avg_n_Wrap(f1->Qle(:,{A:B},{C:D}),0)
Qh_CTL       = dim_avg_n_Wrap(f1->Qh(:,{A:B},{C:D}),0)
Qg_CTL       = dim_avg_n_Wrap(f1->Qg(:,{A:B},{C:D}),0)
SMsurf_CTL   = dim_avg_n_Wrap(f1->SoilMoist(:,0,{A:B},{C:D}),0)
SoilTemp_CTL = dim_avg_n_Wrap(f1->SoilTemp(:,:,{A:B},{C:D}),0)
SoilMoist_CTL= dim_avg_n_Wrap(f1->SoilMoist(:,:,{A:B},{C:D}),0)
LAI_CTL      = dim_avg_n_Wrap(f1->LAI(:,{A:B},{C:D}),0)
GWwb_CTL     = dim_avg_n_Wrap(f1->GWMoist(:,{A:B},{C:D}),0)

Rain_CTL     = 0.
Tair_CTL     = 0.
Qair_CTL     = 0.
Psurf_CTL    = 0.
SWdown_CTL   = 0.
LWdown_CTL   = 0.
Wind_CTL     = 0.
Evap_CTL     = 0.
TVeg_CTL     = 0.
ESoil_CTL    = 0.
ECanop_CTL   = 0.
Runoff_CTL   = 0.
Qrec_CTL     = 0.
Rnet_CTL     = 0.
Qle_CTL      = 0.
Qh_CTL       = 0.
Qg_CTL       = 0.
SoilTemp_CTL = 0.
SoilMoist_CTL= 0.
SMsurf_CTL   = 0.
LAI_CTL      = 0.
GWwb_CTL     = 0.
delete(f1)

filename2 = path2+"/cable_out_"+year_s+".nc"
f2        = addfile (filename2,"r")

Rain_SEN     = dim_sum_n_Wrap(f2->Rainf(:,{A:B},{C:D}),0)
Tair_SEN     = dim_avg_n_Wrap(f2->Tair(:,{A:B},{C:D}),0)
Qair_SEN     = dim_avg_n_Wrap(f2->Qair(:,{A:B},{C:D}),0)
Psurf_SEN    = dim_avg_n_Wrap(f2->PSurf(:,{A:B},{C:D}),0)
SWdown_SEN   = dim_avg_n_Wrap(f2->SWdown(:,{A:B},{C:D}),0)
LWdown_SEN   = dim_avg_n_Wrap(f2->LWdown(:,{A:B},{C:D}),0)
Wind_SEN     = dim_avg_n_Wrap(f2->Wind(:,{A:B},{C:D}),0)

Evap_SEN     = dim_sum_n_Wrap(f2->Evap(:,{A:B},{C:D}),0)
TVeg_SEN     = dim_sum_n_Wrap(f2->TVeg(:,{A:B},{C:D}),0)
ESoil_SEN    = dim_sum_n_Wrap(f2->ESoil(:,{A:B},{C:D}),0)
ECanop_SEN   = dim_sum_n_Wrap(f2->ECanop(:,{A:B},{C:D}),0)
Runoff_SEN   = dim_sum_n_Wrap(f2->Qs(:,{A:B},{C:D}),0)
Qrec_SEN     = dim_sum_n_Wrap(f2->Qrecharge(:,{A:B},{C:D}),0)
Rnet_SEN     = dim_avg_n_Wrap(f2->Rnet(:,{A:B},{C:D}),0)
Qle_SEN      = dim_avg_n_Wrap(f2->Qle(:,{A:B},{C:D}),0)
Qh_SEN       = dim_avg_n_Wrap(f2->Qh(:,{A:B},{C:D}),0)
Qg_SEN       = dim_avg_n_Wrap(f2->Qg(:,{A:B},{C:D}),0)
SMsurf_SEN   = dim_avg_n_Wrap(f2->SoilMoist(:,0,{A:B},{C:D}),0)
SoilTemp_SEN = dim_avg_n_Wrap(f2->SoilTemp(:,:,{A:B},{C:D}),0)
SoilMoist_SEN= dim_avg_n_Wrap(f2->SoilMoist(:,:,{A:B},{C:D}),0)
LAI_SEN      = dim_avg_n_Wrap(f2->LAI(:,{A:B},{C:D}),0)
GWwb_SEN     = dim_avg_n_Wrap(f2->GWMoist(:,{A:B},{C:D}),0)

Rain_SEN     = 0.
Tair_SEN     = 0.
Qair_SEN     = 0.
Psurf_SEN    = 0.
SWdown_SEN   = 0.
LWdown_SEN   = 0.
Wind_SEN     = 0.
Evap_SEN     = 0.
TVeg_SEN     = 0.
ESoil_SEN    = 0.
ECanop_SEN   = 0.
Runoff_SEN   = 0.
Qrec_SEN     = 0.
Rnet_SEN     = 0.
Qle_SEN      = 0.
Qh_SEN       = 0.
Qg_SEN       = 0.
SoilTemp_SEN = 0.
SoilMoist_SEN= 0.
SMsurf_SEN   = 0.
LAI_SEN      = 0.
GWwb_SEN     = 0.
delete(f2)

do year = year_s,year_e
    print("year = " +year)

    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    filename1 = path1+"/cable_out_"+year+".nc"
    f1 = addfile (filename1,"r")

    filename2 = path2+"/cable_out_"+year+".nc"
    f2 = addfile (filename2,"r")

    do month = 1, 12
      Rain_CTL     = Rain_CTL   + f1->Rainf(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      Evap_CTL     = Evap_CTL   + f1->Evap(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      TVeg_CTL     = TVeg_CTL   + f1->TVeg(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      ESoil_CTL    = ESoil_CTL  + f1->ESoil(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      ECanop_CTL   = ECanop_CTL + f1->ECanop(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      Runoff_CTL   = Runoff_CTL + f1->Qs(month-1,{A:B},{C:D})*dom(month-1)*3600*24 \
                                  + f1->Qsb(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      Qrec_CTL     = Qrec_CTL + f1->Qrecharge(month-1,{A:B},{C:D})*dom(month-1)*3600*24

      Rain_SEN     = Rain_SEN   + f2->Rainf(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      Evap_SEN     = Evap_SEN   + f2->Evap(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      TVeg_SEN     = TVeg_SEN   + f2->TVeg(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      ESoil_SEN    = ESoil_SEN  + f2->ESoil(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      ECanop_SEN   = ECanop_SEN + f2->ECanop(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      Runoff_SEN   = Runoff_SEN + f2->Qs(month-1,{A:B},{C:D})*dom(month-1)*3600*24 \
                                  + f2->Qsb(month-1,{A:B},{C:D})*dom(month-1)*3600*24
      Qrec_SEN     = Qrec_SEN + f2->Qrecharge(month-1,{A:B},{C:D})*dom(month-1)*3600*24
    end do ; month

    Tair_CTL     = Tair_CTL   + dim_avg_n_Wrap(f1->Tair(:,{A:B},{C:D}),0)-273.15
    Qair_CTL     = Qair_CTL   + dim_avg_n_Wrap(f1->Qair(:,{A:B},{C:D}),0)
    Psurf_CTL    = Psurf_CTL  + dim_avg_n_Wrap(f1->PSurf(:,{A:B},{C:D}),0)
    SWdown_CTL   = SWdown_CTL + dim_avg_n_Wrap(f1->SWdown(:,{A:B},{C:D}),0)
    LWdown_CTL   = LWdown_CTL + dim_avg_n_Wrap(f1->LWdown(:,{A:B},{C:D}),0)
    Wind_CTL     = Wind_CTL   + dim_avg_n_Wrap(f1->Wind(:,{A:B},{C:D}),0)
    Rnet_CTL     = Rnet_CTL   + dim_avg_n_Wrap(f1->Rnet(:,{A:B},{C:D}),0)
    Qle_CTL      = Qle_CTL    + dim_avg_n_Wrap(f1->Qle(:,{A:B},{C:D}),0)
    Qh_CTL       = Qh_CTL     + dim_avg_n_Wrap(f1->Qh(:,{A:B},{C:D}),0)
    Qg_CTL       = Qg_CTL     + dim_avg_n_Wrap(f1->Qg(:,{A:B},{C:D}),0)
    LAI_CTL      = LAI_CTL    + dim_avg_n_Wrap(f1->LAI(:,{A:B},{C:D}),0)
    SMsurf_CTL   = SMsurf_CTL + ( dim_avg_n_Wrap(f1->SoilMoist(:,0,{A:B},{C:D}),0)*0.022 \
                                + dim_avg_n_Wrap(f1->SoilMoist(:,1,{A:B},{C:D}),0)*0.058 \
                                + dim_avg_n_Wrap(f1->SoilMoist(:,2,{A:B},{C:D}),0)*(0.1 - 0.022 - 0.058)) \
                                / 0.1
    GWwb_CTL      = GWwb_CTL   + dim_avg_n_Wrap(f1->GWMoist(:,{A:B},{C:D}),0)

    Tair_SEN     = Tair_SEN   + dim_avg_n_Wrap(f2->Tair(:,{A:B},{C:D}),0)-273.15
    Qair_SEN     = Qair_SEN   + dim_avg_n_Wrap(f2->Qair(:,{A:B},{C:D}),0)
    Psurf_SEN    = Psurf_SEN  + dim_avg_n_Wrap(f2->PSurf(:,{A:B},{C:D}),0)
    SWdown_SEN   = SWdown_SEN + dim_avg_n_Wrap(f2->SWdown(:,{A:B},{C:D}),0)
    LWdown_SEN   = LWdown_SEN + dim_avg_n_Wrap(f2->LWdown(:,{A:B},{C:D}),0)
    Wind_SEN     = Wind_SEN   + dim_avg_n_Wrap(f2->Wind(:,{A:B},{C:D}),0)
    Rnet_SEN     = Rnet_SEN   + dim_avg_n_Wrap(f2->Rnet(:,{A:B},{C:D}),0)
    Qle_SEN      = Qle_SEN    + dim_avg_n_Wrap(f2->Qle(:,{A:B},{C:D}),0)
    Qh_SEN       = Qh_SEN     + dim_avg_n_Wrap(f2->Qh(:,{A:B},{C:D}),0)
    Qg_SEN       = Qg_SEN     + dim_avg_n_Wrap(f2->Qg(:,{A:B},{C:D}),0)
    LAI_SEN      = LAI_SEN    + dim_avg_n_Wrap(f2->LAI(:,{A:B},{C:D}),0)
    SMsurf_SEN   = SMsurf_SEN + ( dim_avg_n_Wrap(f2->SoilMoist(:,0,{A:B},{C:D}),0)*0.022 \
                                + dim_avg_n_Wrap(f2->SoilMoist(:,1,{A:B},{C:D}),0)*0.058 \
                                + dim_avg_n_Wrap(f2->SoilMoist(:,2,{A:B},{C:D}),0)*(0.1 - 0.022 - 0.058)) \
                                / 0.1
    GWwb_SEN      = GWwb_SEN   + dim_avg_n_Wrap(f2->GWMoist(:,{A:B},{C:D}),0)
    do soil = 0,5
      SoilTemp_CTL(soil,:,:) = SoilTemp_CTL(soil,:,:) + ( dim_avg_n_Wrap(f1->SoilTemp(:,soil,{A:B},{C:D}),0)-273.15)\
                                 *Soil_thickness(soil)/sum(Soil_thickness)
      SoilMoist_CTL(soil,:,:)= SoilMoist_CTL(soil,:,:) + ( dim_avg_n_Wrap(f1->SoilMoist(:,soil,{A:B},{C:D}),0))\
                                 *Soil_thickness(soil)/sum(Soil_thickness)
      SoilTemp_SEN(soil,:,:) = SoilTemp_SEN(soil,:,:) + ( dim_avg_n_Wrap(f2->SoilTemp(:,soil,{A:B},{C:D}),0)-273.15)\
                                 *Soil_thickness(soil)/sum(Soil_thickness)
      SoilMoist_SEN(soil,:,:)= SoilMoist_SEN(soil,:,:) + ( dim_avg_n_Wrap(f2->SoilMoist(:,soil,{A:B},{C:D}),0))\
                                 *Soil_thickness(soil)/sum(Soil_thickness)
    end do ; soil
    delete(f1)
end do ; year
Rain_CTL     = Rain_CTL/year_sum
Evap_CTL     = Evap_CTL/year_sum
TVeg_CTL     = TVeg_CTL/year_sum
ESoil_CTL    = ESoil_CTL/year_sum
ECanop_CTL   = ECanop_CTL/year_sum
Runoff_CTL   = Runoff_CTL/year_sum
Qrec_CTL     = Qrec_CTL/year_sum
Tair_CTL     = Tair_CTL/year_sum
Qair_CTL     = Qair_CTL/year_sum
Psurf_CTL    = Psurf_CTL/year_sum
SWdown_CTL   = SWdown_CTL/year_sum
LWdown_CTL   = LWdown_CTL/year_sum
Wind_CTL     = Wind_CTL/year_sum
Rnet_CTL     = Rnet_CTL/year_sum
Qle_CTL      = Qle_CTL/year_sum
Qh_CTL       = Qh_CTL/year_sum
Qg_CTL       = Qg_CTL/year_sum
LAI_CTL      = LAI_CTL/year_sum
SMsurf_CTL   = SMsurf_CTL/year_sum
SoilTemp_CTL = SoilTemp_CTL/year_sum
SoilMoist_CTL= SoilMoist_CTL/year_sum
GWwb_CTL     = GWwb_CTL/year_sum

Rain_SEN     = Rain_SEN/year_sum
Evap_SEN     = Evap_SEN/year_sum
TVeg_SEN     = TVeg_SEN/year_sum
ESoil_SEN    = ESoil_SEN/year_sum
ECanop_SEN   = ECanop_SEN/year_sum
Runoff_SEN   = Runoff_SEN/year_sum
Qrec_SEN     = Qrec_SEN/year_sum
Tair_SEN     = Tair_SEN/year_sum
Qair_SEN     = Qair_SEN/year_sum
Psurf_SEN    = Psurf_SEN/year_sum
SWdown_SEN   = SWdown_SEN/year_sum
LWdown_SEN   = LWdown_SEN/year_sum
Wind_SEN     = Wind_SEN/year_sum
Rnet_SEN     = Rnet_SEN/year_sum
Qle_SEN      = Qle_SEN/year_sum
Qh_SEN       = Qh_SEN/year_sum
Qg_SEN       = Qg_SEN/year_sum
LAI_SEN      = LAI_SEN/year_sum
SMsurf_SEN   = SMsurf_SEN/year_sum
SoilTemp_SEN = SoilTemp_SEN/year_sum
SoilMoist_SEN= SoilMoist_SEN/year_sum
GWwb_SEN     = GWwb_SEN/year_sum

; ================================ READ GLEAM ==================================
gleam_ET = path_GLEAM+"/E_1980_2018_GLEAM_v3.3a_YR.nc"
gleam_T  = path_GLEAM+"/Et_1980_2018_GLEAM_v3.3a_YR.nc"
gleam_Es = path_GLEAM+"/Eb_1980_2018_GLEAM_v3.3a_YR.nc"
gleam_Ec = path_GLEAM+"/Ei_1980_2018_GLEAM_v3.3a_YR.nc"
gleam_SM = path_GLEAM+"/SMsurf_1980_2018_GLEAM_v3.3a_YR.nc"

G_ET     = addfile (gleam_ET,"r")
G_T      = addfile (gleam_T,"r")
G_Es     = addfile (gleam_Es,"r")
G_Ec     = addfile (gleam_Ec,"r")
G_SM     = addfile (gleam_SM,"r")

E_tmp = G_ET->E(:,{C:D},{A:B})
E_tmp = doubletofloat(E_tmp)
T_tmp = G_T->Et(:,{C:D},{A:B})
T_tmp = doubletofloat(T_tmp)
Es_tmp = G_Es->Eb(:,{C:D},{A:B})
Es_tmp = doubletofloat(Es_tmp)
Ec_tmp = G_Ec->Ei(:,{C:D},{A:B})
Ec_tmp = doubletofloat(Ec_tmp)
SM_tmp = G_SM->SMsurf(:,{C:D},{A:B})
SM_tmp = doubletofloat(SM_tmp)

Evap_GLEAM     = dim_avg_n_Wrap(E_tmp(time|23:28,lat|:, lon|:),0)
TVeg_GLEAM     = dim_avg_n_Wrap(T_tmp(time|23:28,lat|:, lon|:),0)
ESoil_GLEAM    = dim_avg_n_Wrap(Es_tmp(time|23:28,lat|:, lon|:),0)
ECanop_GLEAM   = dim_avg_n_Wrap(Ec_tmp(time|23:28,lat|:, lon|:),0)
SMsurf_GLEAM   = dim_avg_n_Wrap(SM_tmp(time|23:28,lat|:, lon|:),0)
;printVarSummary(Evap_GLEAM)
delete(G_ET)
delete(G_T)
delete(G_Es)
delete(G_Ec)
delete(G_SM)

; ================================ READ CLASS ==================================
file_CLASS = path_CLASS + "/CLASS_v1.1_"+year_s+".nc"
CL = addfile (file_CLASS,"r")

Rnet_CLASS= dim_avg_n_Wrap(CL->rs(:,{A:B},{C:D}),0)
Qle_CLASS = dim_avg_n_Wrap(CL->hfls(:,{A:B},{C:D}),0)
Qh_CLASS  = dim_avg_n_Wrap(CL->hfss(:,{A:B},{C:D}),0)
Qg_CLASS  = dim_avg_n_Wrap(CL->hfds(:,{A:B},{C:D}),0)
delete(CL)
if year_e .gt. 2009 then
  CL_year_end = 2009
else
  CL_year_end = year_e
end if
do year = year_s+1,CL_year_end

    file_CLASS = path_CLASS + "/CLASS_v1.1_"+year+".nc"
    CL = addfile (file_CLASS,"r")

    Rnet_CLASS= Rnet_CLASS + dim_avg_n_Wrap(CL->rs(:,{A:B},{C:D}),0)
    Qle_CLASS = Qle_CLASS + dim_avg_n_Wrap(CL->hfls(:,{A:B},{C:D}),0)
    Qh_CLASS  = Qh_CLASS + dim_avg_n_Wrap(CL->hfss(:,{A:B},{C:D}),0)
    Qg_CLASS  = Qg_CLASS + dim_avg_n_Wrap(CL->hfds(:,{A:B},{C:D}),0)
    delete(CL)
end do ;year
Qle_CLASS  = Qle_CLASS/year_sum
Qh_CLASS   = Qh_CLASS/year_sum
Qg_CLASS   = Qg_CLASS/year_sum


; ================================ READ DOLCE ==================================
file_DOLCE = path_DOLCE +"/DOLCE_v1.0_"+year_s+".nc"
;print(file_DOLCE)
DOLCE      = addfile (file_DOLCE,"r")
Qle_DOLCE  = dim_avg_n_Wrap(DOLCE->ET(:,{A:B},{C:D}),0)

do year = year_s+1,year_e
  file_DOLCE = path_DOLCE +"/DOLCE_v1.0_"+year+".nc"
  DOLCE      = addfile (file_DOLCE,"r")
  Qle_DOLCE  = Qle_DOLCE + dim_avg_n_Wrap(DOLCE->ET(:,{A:B},{C:D}),0)
  delete(DOLCE)
end do ;year
Qle_DOLCE   = Qle_DOLCE/year_sum

; ================================ READ LORA ===================================
file_LORA   = path_LORA +"/LORA_v1.0_"+year_s+".nc"
L           = addfile (file_LORA,"r")
Runoff_LORA = dim_avg_n_Wrap(L->mrro(:,{A:B},{C:D}),0)
Runoff_LORA = 0.

do year = year_s,year_e
  if ( mod(year,4) .eq. 0) then
     dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
  else
     dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  end if

  file_LORA   = path_LORA +"/LORA_v1.0_"+year+".nc"
  L           = addfile (file_LORA,"r")

  do month = 0,11
    Runoff_LORA = Runoff_LORA + L->mrro(month,{A:B},{C:D})*dom(month)*3600*24
  end do
  delete(L)
end do ;year

Runoff_LORA   = Runoff_LORA/year_sum

; ================================ READ CERES ==================================
file_CERES = path_CERES +"/CERES_EBAF-Surface_Ed2.8_Subset_200003-201502.nc"
CERES      = addfile (file_CERES,"r")
Rnet_CERES = dim_avg_n_Wrap(CERES->sfc_net_tot_all_mon(34:105,{A:B},{C:D}),0) ; 2003-2008
; sfc_net_tot_all_mon "Surface Net Total Flux, Monthly Means, All-Sky conditions"
; sfc_net_tot_clr_mon "Surface Net Total Flux, Monthly Means, Clear-Sky conditions"
delete(CERES)


; ================================= Plotting ===================================
res                    = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;��ɫ
res@tmBorderThicknessF = 3.0

res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ����ҳ
res@lbLabelBarOn       = True

;************** ����labelbar ***************
res@lbBoxLinesOn       = True                       ;�ر�lbar box ����
res@lbTitleFont        = 25
res@lbLabelFont        = 25
res@lbTitleFontHeightF = 0.013
res@lbLabelFontHeightF = 0.013
res@txString      = ""
res@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont      = 25

;*************** ���õ�ֵ�� *****************
res@cnLinesOn          = False                       ; �رյ�ֵ������
res@cnLineColor        = "black"
res@cnLineThicknessF   = 1.5
res@cnLineLabelsOn     = False
res@gsnMaximize        = True
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode = "ExplicitLevels"

;************ ����������ֵ��Χ **************
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015
res@tmYLMode  = "Explicit"
res@tmXBMode  = "Explicit"

;**************** ͼ�������� *****************
pres                    = True                                      ; ���� panel resources.
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
;pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01

res@gsnLeftString = ""
res@gsnRightString = ""

; ============================= plot met input =============================
  if (plot_type .eq. "met_force") .or. (plot_type .eq. "all") then

    pic0 = "./plots/LIS_vs_obs_met_"+year_s+"-"+year_e+"_"+case_name
    wks0 = gsn_open_wks("pdf",pic0)
    gsn_define_colormap(wks0,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    ; Rain
    plots = new(3,graphic)
    pres@txString = "Rain (mm/y)"
    levels = fspan(0,1000,21)
    res@cnLevels = levels

    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks0,Rain_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks0,Rain_SEN,res)

    levels = fspan(-100.,100.,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks0,Rain_CTL - Rain_SEN,res)
    gsn_panel(wks0,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; Tair
    pres@txString = "Tair (C)"
    levels = fspan(0,30,16)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks0,Tair_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks0,Tair_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-5.,5.,10)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks0,Tair_CTL - Tair_SEN,res)
    gsn_panel(wks0,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; Qair
    pres@txString = "Qair (%)"
    levels = fspan(0,1,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks0,Qair_CTL*100.,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks0,Qair_SEN*100.,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-0.5,0.5,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks0,Qair_CTL*100. - Qair_SEN*100.,res)
    gsn_panel(wks0,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; Psurf
    pres@txString = "Psurf (hPa)"
    levels = fspan(900,1200,31)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"

    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks0,Psurf_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks0,Psurf_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-100,100,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks0,Psurf_CTL - Psurf_SEN,res)
    gsn_panel(wks0,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; SWdown
    pres@txString = "SWdown (W/m^2)"
    levels = fspan(0,500,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks0,SWdown_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks0,SWdown_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-100,100,41)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks0,SWdown_CTL - SWdown_SEN,res)
    gsn_panel(wks0,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; LWdown
    pres@txString = "LWdown (W/m^2)"
    levels = fspan(0,500,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks0,LWdown_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks0,LWdown_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-100,100,41)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks0,LWdown_CTL - LWdown_SEN,res)
    gsn_panel(wks0,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; Wind
    pres@txString = "Wind (m/s)"
    levels = fspan(0,10,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks0,Wind_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks0,Wind_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-3,3,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks0,Wind_CTL - Wind_SEN,res)
    gsn_panel(wks0,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)
    delete(plots)
  end if


  ; ============================= plot water balance ===========================
  if (plot_type .eq. "water_bal") .or. (plot_type .eq. "all") then
    pic1 = "./plots/LIS_vs_obs_water_balance_"+year_s+"-"+year_e+"_"+case_name
    wks1 = gsn_open_wks("pdf",pic1)
    gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    ; LAI
    plots = new(3,graphic)
    pres@txString = "LAI (mm/y)"
    levels = fspan(0,6,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks1,LAI_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks1,LAI_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-3,3,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks1,LAI_CTL - LAI_SEN,res)
    gsn_panel(wks1,(/plots/),(/1,3/),pres)
    delete(plots)

    ; Evap
    plots = new(4,graphic)
    pres@txString = "Evap (mm/y)"
    levels = fspan(0,1000,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks1,Evap_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks1,Evap_SEN,res)
    res@tiMainString = "GLEAM"
    plots(2) = gsn_csm_contour(wks1,Evap_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-300,300,31)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(3) = gsn_csm_contour(wks1,Evap_CTL - Evap_SEN,res)
    gsn_panel(wks1,(/plots/),(/1,4/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; TVeg
    pres@txString = "TVeg (mm/y)"
    levels = fspan(0,1000,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks1,TVeg_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks1,TVeg_SEN,res)
    res@tiMainString = "GLEAM"
    plots(2) = gsn_csm_contour(wks1,TVeg_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-300,300,31)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(3) = gsn_csm_contour(wks1,TVeg_CTL - TVeg_SEN,res)
    gsn_panel(wks1,(/plots/),(/1,4/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; ESoil
    pres@txString = "ESoil (mm/y)"
    levels = fspan(0,1000,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks1,ESoil_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks1,ESoil_SEN,res)
    res@tiMainString = "GLEAM"
    plots(2) = gsn_csm_contour(wks1,ESoil_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-300,300,31)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(3) = gsn_csm_contour(wks1,ESoil_CTL - ESoil_SEN,res)
    gsn_panel(wks1,(/plots/),(/1,4/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; ECanop
    pres@txString = "ECanop (mm/y)"
    levels = fspan(0,1000,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks1,ECanop_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks1,ECanop_SEN,res)
    res@tiMainString = "GLEAM"
    plots(2) = gsn_csm_contour(wks1,ECanop_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-300,300,31)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(3) = gsn_csm_contour(wks1,ECanop_CTL - ECanop_SEN,res)
    gsn_panel(wks1,(/plots/),(/1,4/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; Runoff
    pres@txString = "Runoff (mm/y)"
    levels = fspan(0,1000,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks1,Runoff_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks1,Runoff_SEN,res)
    res@tiMainString = "LORA"
    plots(2) = gsn_csm_contour(wks1,Runoff_LORA,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-300,300,31)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(3) = gsn_csm_contour(wks1,Runoff_CTL - Runoff_SEN,res)
    gsn_panel(wks1,(/plots/),(/1,4/),pres)
    delete(levels)
    delete(res@cnLevels)
    delete(plots)

    ; Recharge
    plots = new(3,graphic)
    pres@txString = "Recharge (mm/y)"
    levels = fspan(-50,50,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks1,Qrec_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks1,Qrec_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-50,50,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks1,Qrec_CTL - Qrec_SEN,res)
    gsn_panel(wks1,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)
    delete(plots)

    ; SMsurf
    plots = new(4,graphic)
    pres@txString = "SMsurf (m3/m3)"
    levels = fspan(0,0.5,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks1,SMsurf_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks1,SMsurf_SEN,res)
    res@tiMainString = "GLEAM"
    plots(2) = gsn_csm_contour(wks1,SMsurf_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-0.2,0.2,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(3) = gsn_csm_contour(wks1,SMsurf_CTL - SMsurf_SEN,res)
    gsn_panel(wks1,(/plots/),(/1,4/),pres)
    delete(levels)
    delete(res@cnLevels)
    delete(plots)

    ; ; ET/P
    ; plots = new(2,graphic)
    ; pres@txString = "Evap/Rain (-)"
    ; levels = fspan(0,1.,21)
    ; res@cnLevels = levels
    ; res@cnFillPalette = "WhiteBlueGreenYellowRed"
    ; res@tiMainString = "CTL"
    ; plots(0) = gsn_csm_contour(wks1,Evap_CTL/Rain_CTL,res)
    ; res@tiMainString = "SEN"
    ; plots(1) = gsn_csm_contour(wks1,Evap_SEN/Rain_SEN,res)
    ; ; res@tiMainString = "GLEAM"
    ; ; plots(2) = gsn_csm_contour(wks1,Evap_GLEAM/Rain_SEN,res)
    ; ; gsn_panel(wks1,(/plots/),(/1,3/),pres)
    ; gsn_panel(wks1,(/plots/),(/1,2/),pres)
    ;
    ; ; R/ET
    ; pres@txString = "Runoff/Evap (-)"
    ; levels = fspan(0,1.,21)
    ; res@cnLevels = levels
    ; res@cnFillPalette = "WhiteBlueGreenYellowRed"
    ; res@tiMainString = "CTL"
    ; plots(0) = gsn_csm_contour(wks1,Runoff_CTL/Evap_CTL,res)
    ; res@tiMainString = "SEN"
    ; plots(1) = gsn_csm_contour(wks1,Runoff_SEN/Evap_SEN,res)
    ; ; res@tiMainString = "LORA/GLEAM"
    ; ; plots(2) = gsn_csm_contour(wks1,Runoff_LORA/Evap_GLEAM,res)
    ; ; gsn_panel(wks1,(/plots/),(/1,3/),pres)
    ; gsn_panel(wks1,(/plots/),(/1,2/),pres)
    ; delete(plots)
    ; delete(levels)
    ; delete(res@cnLevels)


    plots = new(3,graphic)
    ; Soil Moist
    do soil = 0,5
        pres@txString = "Soil Moisture (m3/m3)"
        levels = fspan(0,0.5,21)
        res@cnLevels = levels
        res@cnFillPalette = "WhiteBlueGreenYellowRed"
        res@tiMainString = "CTL"
        plots(0) = gsn_csm_contour(wks1,SoilMoist_CTL(soil,:,:),res)
        res@tiMainString = "SEN"
        plots(1) = gsn_csm_contour(wks1,SoilMoist_SEN(soil,:,:),res)
        delete(levels)
        delete(res@cnLevels)
        levels = fspan(-0.2,0.2,21)
        res@cnLevels = levels
        res@cnFillPalette = "ViBlGrWhYeOrRe"
        res@tiMainString = "CTL - SEN"
        plots(2) = gsn_csm_contour(wks1,SoilMoist_CTL(soil,:,:) - SoilMoist_SEN(soil,:,:),res)
        gsn_panel(wks1,(/plots/),(/1,3/),pres)
        delete(levels)
        delete(res@cnLevels)
    end do

    pres@txString = "Aquifer Moisture (m3/m3)"
    levels = fspan(0,0.5,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks1,GWwb_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks1,GWwb_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-0.2,0.2,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks1,GWwb_CTL - GWwb_SEN,res)
    gsn_panel(wks1,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)
    delete(plots)

end if

; ============================ plot energy balance =============================
  if (plot_type .eq. "energy_bal") .or. (plot_type .eq. "all") then
    pic2 = "./plots/LIS_vs_obs_energy_balance_"+year_s+"-"+year_e+"_"+case_name
    wks2 = gsn_open_wks("pdf",pic2)
    gsn_define_colormap(wks2,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    ; Rnet
    plots = new(4,graphic)
    pres@txString = "Rnet (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks2,Rnet_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks2,Rnet_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-50,50,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks2,Rnet_CTL - Rnet_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CERES"
    plots(3) = gsn_csm_contour(wks2,Rnet_CERES,res)
    gsn_panel(wks2,(/plots/),(/1,4/),pres)
    delete(plots)

    ; Qle
    plots = new(5,graphic)
    pres@txString = "Qle (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks2,Qle_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks2,Qle_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-50,50,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks2,Qle_CTL - Qle_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@tiMainString = "CLASS"
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    plots(3) = gsn_csm_contour(wks2,Qle_CLASS,res)
    res@tiMainString = "DOLCE"
    plots(4) = gsn_csm_contour(wks2,Qle_DOLCE,res)
    gsn_panel(wks2,(/plots/),(/2,3/),pres)
    delete(plots)

    ; Qh
    plots = new(4,graphic)
    pres@txString = "Qh (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks2,Qh_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks2,Qh_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-50,50,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks2,Qh_CTL - Qh_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CLASS"
    plots(3) = gsn_csm_contour(wks2,Qh_CLASS,res)
    gsn_panel(wks2,(/plots/),(/1,4/),pres)

    ; Qg
    pres@txString = "Qg (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CTL"
    plots(0) = gsn_csm_contour(wks2,Qg_CTL,res)
    res@tiMainString = "SEN"
    plots(1) = gsn_csm_contour(wks2,Qg_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-50,50,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks2,Qg_CTL - Qg_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = "CLASS"
    plots(3) = gsn_csm_contour(wks2,Qg_CLASS,res)
    gsn_panel(wks2,(/plots/),(/1,4/),pres)
  end if

end
