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
    Opt@DstFileName       = "1degree.nc"    ; destination files  ???

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
CTL_name = "ERAI_SE39_CTL"
SEN_name = "ERAI_SE39_Hyds_RUC"

; "ERAI_SE39_CTL"
; "ERAI_SE39_Hyds_RUC"
; "ERAI_SE39_Hydr_param_RUC"
; "ERAI_SE39_Hyre"
; "ERAI_SE39_AQdry"
; "ERAI_SE39_SMdry"
; "ERAI_SE39_SMdry_AQdry"
; "ERAI_SE39_Extrct_5"
; "ERAI_SE39_Extrct_25"
; "ERAI_SE39_Extrct_SMdry_25"
path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+CTL_name+"/LIS_output"
path2 = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+SEN_name+"/LIS_output"

; case_name = ["ERAI_ctl","ERAI_beta_exp","ERAI_gw_off","ERAI_or","ERAI_hyre","ERAI_watmove"]

path_GLEAM= "/g/data/w35/Shared_data/Observations/Global_ET_products/GLEAM_v3_3/3_3a/yearly"
; path_CERES= "/g/data/w35/mm3972/data/CERES" ; CERES_EBAF-Surface_Ed2.8_Subset_200003-201502.nc

; ================================= Predefine ==================================
year_s = 2018
year_e = 2018
year_sum = year_e-year_s+1

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Soil_thickness1 = Soil_thickness
;(/0.005, 0.075, 0.154, 0.409, 1.085, 2.872/) ; LIS-CABLE


; ================================== OPTIONS ===================================
plot_type = "all"
force_name= "ERAI" ;"Princeton";"ERAI"


; ================================ Set Regrid ==================================
interp_method = "neareststod"  ;"neareststod"; "conserve" "patch";"bilinear"

filename  = path1+"/LIS.CABLE."+year_s+"090100.d01.nc"
f = addfile (filename,"r")

lat2d    = f->lat(0,:,:)
lon2d    = f->lon(0,:,:)
lon2d    = where(lon2d.lt.0,360+lon2d,lon2d)

minlon  = ceil(min(lon2d))+0.5 ; 116.5516
minlat  = ceil(min(lat2d))+0.5 ; -44.43771
maxlon  = floor(max(lon2d))-0.5 ; 154.5981
maxlat  = floor(max(lat2d))-0.5 ; -12.5853

; study region
A = minlat
B = maxlat
C = minlon
D = maxlon

; ================================ READ CABLE ==================================
filename = path1+"/LIS.CABLE."+year_s+"010100.d01.nc"
f        = addfile (filename,"r")

Rain_CTL     = dim_sum_n_Wrap(f->Rainf_tavg,0)
Tair_CTL     = dim_avg_n_Wrap(f->Tair_f_inst,0)
Qair_CTL     = dim_avg_n_Wrap(f->Qair_f_inst,0)
Psurf_CTL    = dim_avg_n_Wrap(f->Psurf_f_inst,0)
SWdown_CTL   = dim_avg_n_Wrap(f->SWdown_f_inst,0)
LWdown_CTL   = dim_avg_n_Wrap(f->LWdown_f_inst,0)
Wind_CTL     = dim_avg_n_Wrap(f->Wind_f_inst,0)
Evap_CTL     = dim_sum_n_Wrap(f->Evap_tavg,0)
TVeg_CTL     = dim_sum_n_Wrap(f->TVeg_tavg,0)
ESoil_CTL    = dim_sum_n_Wrap(f->ESoil_tavg,0)
ECanop_CTL   = dim_sum_n_Wrap(f->ECanop_tavg,0)
Runoff_CTL   = dim_sum_n_Wrap(f->Qs_tavg,0)
Rnet_CTL     = dim_avg_n_Wrap(f->Swnet_tavg,0)
Qrec_CTL     = dim_avg_n_Wrap(f->Qrec_tavg,0)
Qle_CTL      = dim_avg_n_Wrap(f->Qle_tavg,0)
Qh_CTL       = dim_avg_n_Wrap(f->Qh_tavg,0)
Qg_CTL       = dim_avg_n_Wrap(f->Qg_tavg,0)
SMsurf_CTL   = dim_avg_n_Wrap(f->SoilMoist_inst(:,0,:,:),0)
SoilTemp_CTL = dim_avg_n_Wrap(f->SoilTemp_inst(:,:,:,:),0)
SoilMoist_CTL= dim_avg_n_Wrap(f->SoilMoist_inst(:,:,:,:),0)
LAI_CTL      = dim_avg_n_Wrap(f->LAI_inst,0)
Gwwb_CTL     = dim_avg_n_Wrap(f->GWwb_inst,0)
Fwsoil_CTL   = dim_avg_n_Wrap(f->fwsoil_tavg,0)


;printVarSummary(SMsurf_CTL)
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
Rnet_CTL     = 0.
Qrec_CTL     = 0.
Qle_CTL      = 0.
Qh_CTL       = 0.
Qg_CTL       = 0.
SoilTemp_CTL = 0.
SoilMoist_CTL= 0.
SMsurf_CTL   = 0.
LAI_CTL      = 0.
Gwwb_CTL     = 0.
Fwsoil_CTL   = 0.
delete(f)


Rain_SEN     = Rain_CTL
Tair_SEN     = Tair_CTL
Qair_SEN     = Qair_CTL
Psurf_SEN    = Psurf_CTL
SWdown_SEN   = SWdown_CTL
LWdown_SEN   = LWdown_CTL
Wind_SEN     = Wind_CTL
Evap_SEN     = Evap_CTL
TVeg_SEN     = TVeg_CTL
ESoil_SEN    = ESoil_CTL
ECanop_SEN   = ECanop_CTL
Runoff_SEN   = Runoff_CTL
Rnet_SEN     = Rnet_CTL
Qrec_SEN     = Qrec_CTL
Qle_SEN      = Qle_CTL
Qh_SEN       = Qh_CTL
Qg_SEN       = Qg_CTL
SMsurf_SEN   = SMsurf_CTL
SoilTemp_SEN = SoilTemp_CTL
SoilMoist_SEN= SoilMoist_CTL
LAI_SEN      = LAI_CTL
Gwwb_SEN     = Gwwb_CTL
Fwsoil_SEN   = Fwsoil_CTL

do year = year_s,year_e
    print("year = " +year)

    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    ; _________________________ read LIS-CABLE ____________________________
    do month = 12, 12
      if month .lt. 10 then
        filename  = path1+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
        filename1 = path2+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
      else
        filename  = path1+"/LIS.CABLE."+year+month+"0100.d01.nc"
        filename1 = path2+"/LIS.CABLE."+year+month+"0100.d01.nc"
      end if
      f  = addfile (filename,"r")
      f1 = addfile (filename1,"r")

      Rain_CTL     = Rain_CTL   + dim_sum_n_Wrap(f->Rainf_tavg,0)*3600*24
      Evap_CTL     = Evap_CTL   + dim_sum_n_Wrap(f->Evap_tavg,0)*3600*24
      TVeg_CTL     = TVeg_CTL   + dim_sum_n_Wrap(f->TVeg_tavg,0)*3600*24
      ESoil_CTL    = ESoil_CTL  + dim_sum_n_Wrap(f->ESoil_tavg,0)*3600*24
      ECanop_CTL   = ECanop_CTL + dim_sum_n_Wrap(f->ECanop_tavg,0)*3600*24
      Runoff_CTL   = Runoff_CTL + dim_sum_n_Wrap(f->Qs_tavg,0)*3600*24 \
                                + dim_sum_n_Wrap(f->Qsb_tavg,0)*3600*24
      Qrec_CTL     = Qrec_CTL   + dim_sum_n_Wrap(f->Qrec_tavg,0)*3600*24
      Tair_CTL     = Tair_CTL   + dim_avg_n_Wrap(f->Tair_f_inst,0)-273.15
      Qair_CTL     = Qair_CTL   + dim_avg_n_Wrap(f->Qair_f_inst,0)
      Psurf_CTL    = Psurf_CTL  + dim_avg_n_Wrap(f->Psurf_f_inst,0)/100.
      SWdown_CTL   = SWdown_CTL + dim_avg_n_Wrap(f->SWdown_f_inst,0)
      LWdown_CTL   = LWdown_CTL + dim_avg_n_Wrap(f->LWdown_f_inst,0)
      Wind_CTL     = Wind_CTL   + dim_avg_n_Wrap(f->Wind_f_inst,0)
      Rnet_CTL     = Rnet_CTL   + dim_avg_n_Wrap(f->Swnet_tavg,0)\
                                + dim_avg_n_Wrap(f->Lwnet_tavg,0)
      Qle_CTL      = Qle_CTL    + dim_avg_n_Wrap(f->Qle_tavg,0)
      Qh_CTL       = Qh_CTL     + dim_avg_n_Wrap(f->Qh_tavg,0)
      Qg_CTL       = Qg_CTL     + dim_avg_n_Wrap(f->Qg_tavg,0)
      LAI_CTL      = LAI_CTL    + dim_avg_n_Wrap(f->LAI_inst,0)
      SMsurf_CTL   = SMsurf_CTL + ( dim_avg_n_Wrap(f->SoilMoist_inst(:,0,:,:),0)*0.005 \
                                  + dim_avg_n_Wrap(f->SoilMoist_inst(:,1,:,:),0)*0.075 \
                                  + dim_avg_n_Wrap(f->SoilMoist_inst(:,2,:,:),0)*(0.1 - 0.005 - 0.075)) \
                                  / 0.1
      Gwwb_CTL     = Gwwb_CTL    + dim_avg_n_Wrap(f->GWwb_inst,0);*25.*1000.
      Fwsoil_CTL   = Fwsoil_CTL + dim_avg_n_Wrap(f->fwsoil_tavg,0)

      Rain_SEN     = Rain_SEN   + dim_sum_n_Wrap(f1->Rainf_tavg,0)*3600*24
      Evap_SEN     = Evap_SEN   + dim_sum_n_Wrap(f1->Evap_tavg,0)*3600*24
      TVeg_SEN     = TVeg_SEN   + dim_sum_n_Wrap(f1->TVeg_tavg,0)*3600*24
      ESoil_SEN    = ESoil_SEN  + dim_sum_n_Wrap(f1->ESoil_tavg,0)*3600*24
      ECanop_SEN   = ECanop_SEN + dim_sum_n_Wrap(f1->ECanop_tavg,0)*3600*24
      Runoff_SEN   = Runoff_SEN + dim_sum_n_Wrap(f1->Qs_tavg,0)*3600*24 \
                                + dim_sum_n_Wrap(f1->Qsb_tavg,0)*3600*24
      Qrec_SEN     = Qrec_SEN   + dim_sum_n_Wrap(f1->Qrec_tavg,0)*3600*24
      Tair_SEN     = Tair_SEN   + dim_avg_n_Wrap(f1->Tair_f_inst,0)-273.15
      Qair_SEN     = Qair_SEN   + dim_avg_n_Wrap(f1->Qair_f_inst,0)
      Psurf_SEN    = Psurf_SEN  + dim_avg_n_Wrap(f1->Psurf_f_inst,0)/100.
      SWdown_SEN   = SWdown_SEN + dim_avg_n_Wrap(f1->SWdown_f_inst,0)
      LWdown_SEN   = LWdown_SEN + dim_avg_n_Wrap(f1->LWdown_f_inst,0)
      Wind_SEN     = Wind_SEN   + dim_avg_n_Wrap(f1->Wind_f_inst,0)
      Rnet_SEN     = Rnet_SEN   + dim_avg_n_Wrap(f1->Swnet_tavg,0)\
                                + dim_avg_n_Wrap(f1->Lwnet_tavg,0)
      Qle_SEN      = Qle_SEN    + dim_avg_n_Wrap(f1->Qle_tavg,0)
      Qh_SEN       = Qh_SEN     + dim_avg_n_Wrap(f1->Qh_tavg,0)
      Qg_SEN       = Qg_SEN     + dim_avg_n_Wrap(f1->Qg_tavg,0)
      LAI_SEN      = LAI_SEN    + dim_avg_n_Wrap(f1->LAI_inst,0)
      SMsurf_SEN   = SMsurf_SEN + ( dim_avg_n_Wrap(f1->SoilMoist_inst(:,0,:,:),0)*0.005 \
                                  + dim_avg_n_Wrap(f1->SoilMoist_inst(:,1,:,:),0)*0.075 \
                                  + dim_avg_n_Wrap(f1->SoilMoist_inst(:,2,:,:),0)*(0.1 - 0.005 - 0.075)) \
                                  / 0.1
      Gwwb_SEN     = Gwwb_SEN    + dim_avg_n_Wrap(f1->GWwb_inst,0);*25.*1000.
      Fwsoil_SEN   = Fwsoil_SEN + dim_avg_n_Wrap(f1->fwsoil_tavg,0)

      do soil = 0,5
        SoilTemp_CTL(soil,:,:) = SoilTemp_CTL(soil,:,:) + (dim_avg_n_Wrap(f->SoilTemp_inst(:,soil,:,:),0)-273.15)
                                 ;*Soil_thickness1(soil)/sum(Soil_thickness1)
        SoilMoist_CTL(soil,:,:)= SoilMoist_CTL(soil,:,:) +dim_avg_n_Wrap(f->SoilMoist_inst(:,soil,:,:),0);\
                                 ;*Soil_thickness1(soil)*1000.
        SoilTemp_SEN(soil,:,:) = SoilTemp_SEN(soil,:,:) + (dim_avg_n_Wrap(f1->SoilTemp_inst(:,soil,:,:),0)-273.15)
                                  ;  *Soil_thickness1(soil)/sum(Soil_thickness1)
        SoilMoist_SEN(soil,:,:)= SoilMoist_SEN(soil,:,:) +dim_avg_n_Wrap(f1->SoilMoist_inst(:,soil,:,:),0);\
                                 ;*Soil_thickness1(soil)*1000.
      end do ; soil
      delete(f)
      delete(f1)
    end do ; month
end do ; year

Rain_CTL     = Rain_CTL/year_sum
Evap_CTL     = Evap_CTL/year_sum
TVeg_CTL     = TVeg_CTL/year_sum
ESoil_CTL    = ESoil_CTL/year_sum
ECanop_CTL   = ECanop_CTL/year_sum
Runoff_CTL   = Runoff_CTL/year_sum
Qrec_CTL     = Qrec_CTL/year_sum
Tair_CTL     = Tair_CTL/year_sum;/12
Qair_CTL     = Qair_CTL/year_sum;/12
Psurf_CTL    = Psurf_CTL/year_sum;/12
SWdown_CTL   = SWdown_CTL/year_sum;/12
LWdown_CTL   = LWdown_CTL/year_sum;/12
Wind_CTL     = Wind_CTL/year_sum;/12
Rnet_CTL     = Rnet_CTL/year_sum;/12
Qle_CTL      = Qle_CTL/year_sum;/12
Qh_CTL       = Qh_CTL/year_sum;/12
Qg_CTL       = Qg_CTL/year_sum;/12
SMsurf_CTL   = SMsurf_CTL/year_sum;/12
LAI_CTL      = LAI_CTL/year_sum;/12
SoilTemp_CTL = SoilTemp_CTL/year_sum;/12
SoilMoist_CTL= SoilMoist_CTL/year_sum;/12
Gwwb_CTL     = Gwwb_CTL/year_sum;/12
Fwsoil_CTL   = Fwsoil_CTL/year_sum;/12

Rain_SEN     = Rain_SEN/year_sum
Evap_SEN     = Evap_SEN/year_sum
TVeg_SEN     = TVeg_SEN/year_sum
ESoil_SEN    = ESoil_SEN/year_sum
ECanop_SEN   = ECanop_SEN/year_sum
Runoff_SEN   = Runoff_SEN/year_sum
Qrec_SEN     = Qrec_SEN/year_sum
Tair_SEN     = Tair_SEN/year_sum;/12
Qair_SEN     = Qair_SEN/year_sum;/12
Psurf_SEN    = Psurf_SEN/year_sum;/12
SWdown_SEN   = SWdown_SEN/year_sum;/12
LWdown_SEN   = LWdown_SEN/year_sum;/12
Wind_SEN     = Wind_SEN/year_sum;/12
Rnet_SEN     = Rnet_SEN/year_sum;/12
Qle_SEN      = Qle_SEN/year_sum;/12
Qh_SEN       = Qh_SEN/year_sum;/12
Qg_SEN       = Qg_SEN/year_sum;/12
SMsurf_SEN   = SMsurf_SEN/year_sum;/12
LAI_SEN      = LAI_SEN/year_sum;/12
SoilTemp_SEN = SoilTemp_SEN/year_sum;/12
SoilMoist_SEN= SoilMoist_SEN/year_sum;/12
Gwwb_SEN     = Gwwb_SEN/year_sum;/12
Fwsoil_SEN   = Fwsoil_SEN/year_sum;/12

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

Evap_GLEAM     = dim_avg_n_Wrap(E_tmp(time|:,lat|:, lon|:),0)
TVeg_GLEAM     = dim_avg_n_Wrap(T_tmp(time|:,lat|:, lon|:),0)
ESoil_GLEAM    = dim_avg_n_Wrap(Es_tmp(time|:,lat|:, lon|:),0)
ECanop_GLEAM   = dim_avg_n_Wrap(Ec_tmp(time|:,lat|:, lon|:),0)
SMsurf_GLEAM   = dim_avg_n_Wrap(SM_tmp(time|:,lat|:, lon|:),0)
;printVarSummary(Evap_GLEAM)
delete(G_ET)
delete(G_T)
delete(G_Es)
delete(G_Ec)
delete(G_SM)

; ; ================================ READ CERES ==================================
; file_CERES = path_CERES +"/CERES_EBAF-Surface_Ed2.8_Subset_200003-201502.nc"
; CERES      = addfile (file_CERES,"r")
; Rnet_CERES = dim_avg_n_Wrap(CERES->sfc_net_tot_all_mon(34:105,{A:B},{C:D}),0) ; 2003-2008
; ; sfc_net_tot_all_mon "Surface Net Total Flux, Monthly Means, All-Sky conditions"
; ; sfc_net_tot_clr_mon "Surface Net Total Flux, Monthly Means, Clear-Sky conditions"
; delete(CERES)


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

    pic0 = "./plots/LIS_met_force_"+force_name+"_"+year_s+"-"+year_e
    wks0 = gsn_open_wks("pdf",pic0)
    gsn_define_colormap(wks0,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    plots = new(7,graphic)
    pres@txString = "met force"

    ; Rain
    res@tiMainString = "Rain (mm/y)"
    levels = fspan(0,1000,21)
    res@cnLevels = levels
    plots(0) = gsn_csm_contour(wks0,Rain_CTL,res)
    delete(levels)
    delete(res@cnLevels)

    ; Tair
    res@tiMainString = "Tair (C)"
    levels = fspan(0,30,16)
    res@cnLevels = levels
    plots(1) = gsn_csm_contour(wks0,Tair_CTL,res)
    delete(levels)
    delete(res@cnLevels)

    ; Qair
    res@tiMainString = "Qair (%)"
    levels = fspan(0,1,21)
    res@cnLevels = levels
    plots(2) = gsn_csm_contour(wks0,Qair_CTL*100.,res)
    delete(levels)
    delete(res@cnLevels)

    ; Psurf
    res@tiMainString = "Psurf (hPa)"
    levels = fspan(900,1200,31)
    res@cnLevels = levels
    plots(3) = gsn_csm_contour(wks0,Psurf_CTL,res)
    delete(levels)
    delete(res@cnLevels)

    ; SWdown
    res@tiMainString = "SWdown (W/m^2)"
    levels = fspan(0,500,51)
    res@cnLevels = levels
    plots(4) = gsn_csm_contour(wks0,SWdown_CTL,res)
    delete(levels)
    delete(res@cnLevels)

    ; LWdown
    res@tiMainString = "LWdown (W/m^2)"
    levels = fspan(0,500,51)
    res@cnLevels = levels
    plots(5) = gsn_csm_contour(wks0,LWdown_CTL,res)
    delete(levels)
    delete(res@cnLevels)

    ; Wind
    res@tiMainString = "Wind (m/s)"
    levels = fspan(0,10,21)
    res@cnLevels = levels
    plots(6) = gsn_csm_contour(wks0,Wind_CTL,res)
    gsn_panel(wks0,(/plots/),(/3,3/),pres)
    delete(levels)
    delete(res@cnLevels)
    delete(plots)
  end if


  ; ============================= plot water balance ===========================
  if (plot_type .eq. "water_bal") .or. (plot_type .eq. "all") then
    pic1 = "./plots/LIS_water_balance_"+CTL_name+"_vs_"+SEN_name+"_"+force_name+"_"+year_s+"-"+year_e
    wks1 = gsn_open_wks("pdf",pic1)
    gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    ; LAI
    plots = new(3,graphic)
    pres@txString = "LAI (mm/y)"
    levels = fspan(0,6,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks1,LAI_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks1,LAI_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    levels = fspan(-3,3,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = "CTL - SEN"
    plots(2) = gsn_csm_contour(wks1,LAI_SEN - LAI_CTL,res)
    gsn_panel(wks1,(/plots/),(/1,3/),pres)
    delete(plots)

    ; Evap
    plots = new(4,graphic)
    pres@txString = "Evap (mm/y)"
    ;levels = fspan(0,1000,21)
    levels = fspan(0,500,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks1,Evap_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks1,Evap_SEN,res)
    res@tiMainString = "GLEAM"
    plots(2) = gsn_csm_contour(wks1,Evap_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    ;levels = fspan(-50,50,21)
    levels = fspan(-20,20,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(3) = gsn_csm_contour(wks1,Evap_SEN - Evap_CTL,res)
    gsn_panel(wks1,(/plots/),(/1,4/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; TVeg
    pres@txString = "TVeg (mm/y)"
    ;levels = fspan(0,1000,21)
    levels = fspan(0,500,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks1,TVeg_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks1,TVeg_SEN,res)
    res@tiMainString = "GLEAM"
    plots(2) = gsn_csm_contour(wks1,TVeg_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    ;levels = fspan(-50,50,21)
    levels = fspan(-20,20,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(3) = gsn_csm_contour(wks1,TVeg_SEN - TVeg_CTL,res)
    gsn_panel(wks1,(/plots/),(/1,4/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; ESoil
    pres@txString = "ESoil (mm/y)"
    ;levels = fspan(0,1000,21)
    levels = fspan(0,500,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks1,ESoil_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks1,ESoil_SEN,res)
    res@tiMainString = "GLEAM"
    plots(2) = gsn_csm_contour(wks1,ESoil_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    ;levels = fspan(-50,50,21)
    levels = fspan(-20,20,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(3) = gsn_csm_contour(wks1,ESoil_SEN - ESoil_CTL,res)
    gsn_panel(wks1,(/plots/),(/1,4/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; ECanop
    pres@txString = "ECanop (mm/y)"
    ;levels = fspan(0,1000,21)
    levels = fspan(0,500,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks1,ECanop_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks1,ECanop_SEN,res)
    res@tiMainString = "GLEAM"
    plots(2) = gsn_csm_contour(wks1,ECanop_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    ;levels = fspan(-50,50,21)
    levels = fspan(-20,20,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(3) = gsn_csm_contour(wks1,ECanop_SEN - ECanop_CTL,res)
    gsn_panel(wks1,(/plots/),(/1,4/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; SMsurf
    pres@txString = "SMsurf (m3/m3)"
    ;levels = fspan(0,0.5,21)
    levels = fspan(0,0.2,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks1,SMsurf_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks1,SMsurf_SEN,res)
    res@tiMainString = "GLEAM"
    plots(2) = gsn_csm_contour(wks1,SMsurf_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    ;levels = fspan(-0.1,0.1,21)
    levels = fspan(-0.05,0.05,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(3) = gsn_csm_contour(wks1,SMsurf_SEN - SMsurf_CTL,res)
    gsn_panel(wks1,(/plots/),(/1,4/),pres)
    delete(levels)
    delete(res@cnLevels)
    delete(plots)

    ; Runoff
    plots = new(3,graphic)
    pres@txString = "Runoff (mm/y)"
    ;levels = fspan(0,1000,21)
    levels = fspan(0,500,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks1,Runoff_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks1,Runoff_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    ;levels = fspan(-300,300,31)
    levels = fspan(-50,50,31)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(2) = gsn_csm_contour(wks1,Runoff_SEN - Runoff_CTL,res)
    gsn_panel(wks1,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)
    delete(plots)

    ; Qrec
    plots = new(3,graphic)
    pres@txString = "Recharge (mm/y)"
    levels = fspan(-50,50,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks1,Qrec_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks1,Qrec_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    ;levels = fspan(-20,20,21)
    levels = fspan(-10,10,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(2) = gsn_csm_contour(wks1,Qrec_SEN - Qrec_CTL,res)
    gsn_panel(wks1,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)
    delete(plots)

    ; Fwsoil
    plots = new(3,graphic)
    pres@txString = "Fwsoil (-)"
    levels = fspan(0,1,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks1,Fwsoil_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks1,Fwsoil_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    ;levels = fspan(-20,20,21)
    levels = fspan(-0.2,0.2,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(2) = gsn_csm_contour(wks1,Fwsoil_SEN - Fwsoil_CTL,res)
    gsn_panel(wks1,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)
    delete(plots)


    ; ET/P
    plots = new(2,graphic)
    pres@txString = "Evap/Rain (-)"
    levels = fspan(0,1.,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks1,Evap_CTL/Rain_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks1,Evap_SEN/Rain_SEN,res)
    gsn_panel(wks1,(/plots/),(/1,2/),pres)

    ; ; R/ET
    ; pres@txString = "Runoff/Evap (-)"
    ; levels = fspan(0,1.,21)
    ; res@cnLevels = levels
    ; res@cnFillPalette = "WhiteBlueGreenYellowRed"
    ; res@tiMainString = CTL_name
    ; plots(0) = gsn_csm_contour(wks1,Runoff_CTL/Evap_CTL,res)
    ; res@tiMainString = SEN_name
    ; plots(1) = gsn_csm_contour(wks1,Runoff_SEN/Evap_SEN,res)
    ; gsn_panel(wks1,(/plots/),(/1,2/),pres)
    delete(plots)
    ; delete(levels)
    ; delete(res@cnLevels)


    ; Soil Moist
    do layer=0,5
        plots = new(3,graphic)
        pres@txString = "Soil Water Storage at "+ layer +"-layer (mm)"
        levels = fspan(0,0.5,21)
        res@cnLevels = levels
        res@cnFillPalette = "WhiteBlueGreenYellowRed"
        res@tiMainString = CTL_name
        plots(0) = gsn_csm_contour(wks1,SoilMoist_CTL(layer,:,:),res)
        res@tiMainString = SEN_name
        plots(1) = gsn_csm_contour(wks1,SoilMoist_SEN(layer,:,:),res)
        delete(levels)
        delete(res@cnLevels)
        ;levels = fspan(-0.2,0.2,21)
        levels = fspan(-0.1,0.1,21)
        res@cnLevels = levels
        res@cnFillPalette = "ViBlGrWhYeOrRe"
        res@tiMainString = SEN_name+" - "+ CTL_name
        plots(2) = gsn_csm_contour(wks1,SoilMoist_SEN(layer,:,:) - SoilMoist_CTL(layer,:,:),res)
        gsn_panel(wks1,(/plots/),(/1,3/),pres)
        delete(levels)
        delete(res@cnLevels)
        delete(plots)
    end do

    plots = new(3,graphic)
    pres@txString = "Aquifer Water storage (mm)"
    levels = fspan(0,0.5,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks1,Gwwb_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks1,Gwwb_SEN,res)
    delete(levels)
    delete(res@cnLevels)
    ;levels = fspan(-0.2,0.2,21)
    levels = fspan(-0.1,0.1,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(2) = gsn_csm_contour(wks1,Gwwb_SEN - Gwwb_CTL,res)
    gsn_panel(wks1,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)
    delete(plots)
end if

; ============================ plot energy balance =============================
  if (plot_type .eq. "energy_bal") .or. (plot_type .eq. "all") then
    pic2 = "./plots/LIS_energy_balance_"+CTL_name+"_vs_"+SEN_name+"_"+force_name+"_"+year_s+"-"+year_e
    wks2 = gsn_open_wks("pdf",pic2)
    gsn_define_colormap(wks2,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    ; Rnet
    plots = new(3,graphic)
    pres@txString = "Rnet (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks2,Rnet_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks2,Rnet_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-20,20,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(2) = gsn_csm_contour(wks2,Rnet_SEN - Rnet_CTL,res)
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks2,(/plots/),(/1,3/),pres)

    ; Qle
    pres@txString = "Qle (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks2,Qle_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks2,Qle_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-20,20,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(2) = gsn_csm_contour(wks2,Qle_SEN - Qle_CTL,res)
    gsn_panel(wks2,(/plots/),(/1,3/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; Qh
    pres@txString = "Qh (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks2,Qh_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks2,Qh_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-20,20,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(2) = gsn_csm_contour(wks2,Qh_SEN - Qh_CTL,res)
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks2,(/plots/),(/1,3/),pres)

    ; Qg
    pres@txString = "Qg (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    res@tiMainString = CTL_name
    plots(0) = gsn_csm_contour(wks2,Qg_CTL,res)
    res@tiMainString = SEN_name
    plots(1) = gsn_csm_contour(wks2,Qg_SEN,res)
    delete(levels)
    delete(res@cnLevels)

    levels = fspan(-20,20,21)
    res@cnLevels = levels
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    res@tiMainString = SEN_name+" - "+ CTL_name
    plots(2) = gsn_csm_contour(wks2,Qg_SEN - Qg_CTL,res)
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks2,(/plots/),(/1,3/),pres)
  end if

end
