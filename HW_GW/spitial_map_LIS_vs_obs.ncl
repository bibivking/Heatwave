;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

; =================================== Path =====================================
path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/nuwrf_25km/LISWRF_configs/LIS_output"
path2 = "/g/data/w35/mm3972/model/cable/runs/my_version/outputs"

path_GLEAM= "/g/data/w35/Shared_data/Observations/Global_ET_products/GLEAM_v3_3/3_3a/yearly";
path_CLASS= "/g/data/w35/mm3972/data/CLASS_v1"
path_DOLCE= "/g/data/w35/mm3972/data/DOLCE" ; DOLCE_v1.0_2000.nc
path_LORA = "/g/data/w35/mm3972/data/LORA"  ; LORA_v1.0_2000.nc
path_CERES= "/g/data/w35/mm3972/data/CERES" ; CERES_EBAF-Surface_Ed2.8_Subset_200003-201502.nc

; ================================= Predefine ==================================
year_s = 2003
year_e = 2008
year_sum = year_e-year_s+1

; study region
A = -43.4377
B = -13.5853
C = 117.552
D = 153.598

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

; ================================ READ CABLE ==================================
filename = path1+"/LIS.CABLE."+year_s+"010100.d01.nc"
f        = addfile (filename,"r")

Rain_LIS     = dim_sum_n_Wrap(f->Rainf_tavg,0)
Evap_LIS     = dim_sum_n_Wrap(f->Evap_tavg,0)
TVeg_LIS     = dim_sum_n_Wrap(f->TVeg_tavg,0)
ESoil_LIS    = dim_sum_n_Wrap(f->ESoil_tavg,0)
ECanop_LIS   = dim_sum_n_Wrap(f->ECanop_tavg,0)
Runoff_LIS   = dim_sum_n_Wrap(f->Qs_tavg,0)
Rnet_LIS     = dim_avg_n_Wrap(f->Swnet_tavg,0)
Qle_LIS      = dim_avg_n_Wrap(f->Qle_tavg,0)
Qh_LIS       = dim_avg_n_Wrap(f->Qh_tavg,0)
Qg_LIS       = dim_avg_n_Wrap(f->Qg_tavg,0)
SMsurf_LIS   = dim_avg_n_Wrap(f->SoilMoist_inst(:,0,:,:),0)
SoilTemp_LIS = dim_avg_n_Wrap(f->SoilTemp_inst(:,0,:,:),0)
SoilMoist_LIS= dim_avg_n_Wrap(f->SoilMoist_inst(:,0,:,:),0)
LAI_LIS      = dim_avg_n_Wrap(f->LAI_inst,0)


;printVarSummary(SMsurf_LIS)
Rain_LIS     = 0.
Evap_LIS     = 0.
TVeg_LIS     = 0.
ESoil_LIS    = 0.
ECanop_LIS   = 0.
Runoff_LIS   = 0.
Rnet_LIS     = 0.
Qle_LIS      = 0.
Qh_LIS       = 0.
Qg_LIS       = 0.
SoilTemp_LIS = 0.
SoilMoist_LIS= 0.
SMsurf_LIS   = 0.
LAI_LIS      = 0.
delete(f)

filename1 = path2+"/cable_out_"+year_s+".nc"
f1        = addfile (filename1,"r")

Rain_HESS     = dim_sum_n_Wrap(f1->Rainf(:,{A:B},{C:D}),0)
Evap_HESS     = dim_sum_n_Wrap(f1->Evap(:,{A:B},{C:D}),0)
TVeg_HESS     = dim_sum_n_Wrap(f1->TVeg(:,{A:B},{C:D}),0)
ESoil_HESS    = dim_sum_n_Wrap(f1->ESoil(:,{A:B},{C:D}),0)
ECanop_HESS   = dim_sum_n_Wrap(f1->ECanop(:,{A:B},{C:D}),0)
Runoff_HESS   = dim_sum_n_Wrap(f1->Qs(:,{A:B},{C:D}),0)
Rnet_HESS     = dim_avg_n_Wrap(f1->Rnet(:,{A:B},{C:D}),0)
Qle_HESS      = dim_avg_n_Wrap(f1->Qle(:,{A:B},{C:D}),0)
Qh_HESS       = dim_avg_n_Wrap(f1->Qh(:,{A:B},{C:D}),0)
Qg_HESS       = dim_avg_n_Wrap(f1->Qg(:,{A:B},{C:D}),0)
SMsurf_HESS   = dim_avg_n_Wrap(f1->SoilMoist(:,0,{A:B},{C:D}),0)
SoilTemp_HESS = dim_avg_n_Wrap(f1->SoilTemp(:,0,{A:B},{C:D}),0)
SoilMoist_HESS= dim_avg_n_Wrap(f1->SoilMoist(:,0,{A:B},{C:D}),0)
LAI_HESS      = dim_avg_n_Wrap(f1->LAI(:,{A:B},{C:D}),0)

Rain_HESS     = 0.
Evap_HESS     = 0.
TVeg_HESS     = 0.
ESoil_HESS    = 0.
ECanop_HESS   = 0.
Runoff_HESS   = 0.
Rnet_HESS     = 0.
Qle_HESS      = 0.
Qh_HESS       = 0.
Qg_HESS       = 0.
SoilTemp_HESS = 0.
SoilMoist_HESS= 0.
SMsurf_HESS   = 0.
LAI_HESS      = 0.
delete(f1)

do year = year_s,year_e
    print("year = " +year)
    ; _________________________ read LIS-CABLE ____________________________
    do month = 1, 12
      if month .lt. 10 then
        filename = path1+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
      else
        filename = path1+"/LIS.CABLE."+year+month+"0100.d01.nc"
      end if
      f = addfile (filename,"r")

      Rain_LIS     = Rain_LIS   + dim_sum_n_Wrap(f->Rainf_tavg,0)*3600*24
      Evap_LIS     = Evap_LIS   + dim_sum_n_Wrap(f->Evap_tavg,0)*3600*24
      TVeg_LIS     = TVeg_LIS   + dim_sum_n_Wrap(f->TVeg_tavg,0)*3600*24
      ESoil_LIS    = ESoil_LIS  + dim_sum_n_Wrap(f->ESoil_tavg,0)*3600*24
      ECanop_LIS   = ECanop_LIS + dim_sum_n_Wrap(f->ECanop_tavg,0)*3600*24
      Runoff_LIS   = Runoff_LIS + dim_sum_n_Wrap(f->Qs_tavg,0)*3600*24 \
                                + dim_sum_n_Wrap(f->Qsb_tavg,0)*3600*24
      Rnet_LIS     = Rnet_LIS   + dim_avg_n_Wrap(f->Swnet_tavg,0)\
                                + dim_avg_n_Wrap(f->Lwnet_tavg,0)
      Qle_LIS      = Qle_LIS    + dim_avg_n_Wrap(f->Qle_tavg,0)
      Qh_LIS       = Qh_LIS     + dim_avg_n_Wrap(f->Qh_tavg,0)
      Qg_LIS       = Qg_LIS     + dim_avg_n_Wrap(f->Qg_tavg,0)
      LAI_LIS      = LAI_LIS    + dim_avg_n_Wrap(f->LAI_inst,0)
      SMsurf_LIS   = SMsurf_LIS + ( dim_avg_n_Wrap(f->SoilMoist_inst(:,0,:,:),0)*0.022 \
                                  + dim_avg_n_Wrap(f->SoilMoist_inst(:,1,:,:),0)*0.058 \
                                  + dim_avg_n_Wrap(f->SoilMoist_inst(:,2,:,:),0)*(0.1 - 0.022 - 0.058)) \
                                  / 0.1
      do soil = 0,5
        SoilTemp_LIS = SoilTemp_LIS + (dim_avg_n_Wrap(f->SoilTemp_inst(:,soil,:,:),0)-273.15)\
                                    * Soil_thickness(soil)/sum(Soil_thickness)
        SoilMoist_LIS= SoilMoist_LIS +dim_avg_n_Wrap(f->SoilMoist_inst(:,soil,:,:),0)\
                                       *Soil_thickness(soil)/sum(Soil_thickness)
      end do ; soil
      delete(f)
    end do ; month

    ; _________________________ read CABLE (HESS) _________________________
    filename1 = path2+"/cable_out_"+year+".nc"
    f1 = addfile (filename1,"r")
    Rain_HESS     = Rain_HESS   + dim_sum_n_Wrap(f1->Rainf(:,{A:B},{C:D}),0)*3600*24
    Evap_HESS     = Evap_HESS   + dim_sum_n_Wrap(f1->Evap(:,{A:B},{C:D}),0)*3600*24
    TVeg_HESS     = TVeg_HESS   + dim_sum_n_Wrap(f1->TVeg(:,{A:B},{C:D}),0)*3600*24
    ESoil_HESS    = ESoil_HESS  + dim_sum_n_Wrap(f1->ESoil(:,{A:B},{C:D}),0)*3600*24
    ECanop_HESS   = ECanop_HESS + dim_sum_n_Wrap(f1->ECanop(:,{A:B},{C:D}),0)*3600*24
    Runoff_HESS   = Runoff_HESS + dim_sum_n_Wrap(f1->Qs(:,{A:B},{C:D}),0)*3600*24 \
                                + dim_sum_n_Wrap(f1->Qsb(:,{A:B},{C:D}),0)*3600*24
    Rnet_HESS     = Rnet_HESS   + dim_avg_n_Wrap(f1->Rnet(:,{A:B},{C:D}),0)
    Qle_HESS      = Qle_HESS    + dim_avg_n_Wrap(f1->Qle(:,{A:B},{C:D}),0)
    Qh_HESS       = Qh_HESS     + dim_avg_n_Wrap(f1->Qh(:,{A:B},{C:D}),0)
    Qg_HESS       = Qg_HESS     + dim_avg_n_Wrap(f1->Qg(:,{A:B},{C:D}),0)
    LAI_HESS      = LAI_HESS    + dim_avg_n_Wrap(f1->LAI(:,{A:B},{C:D}),0)
    SMsurf_HESS   = SMsurf_HESS + ( dim_avg_n_Wrap(f1->SoilMoist(:,0,{A:B},{C:D}),0)*0.022 \
                                + dim_avg_n_Wrap(f1->SoilMoist(:,1,{A:B},{C:D}),0)*0.058 \
                                + dim_avg_n_Wrap(f1->SoilMoist(:,2,{A:B},{C:D}),0)*(0.1 - 0.022 - 0.058)) \
                                / 0.1
    do soil = 0,5
      SoilTemp_HESS = SoilTemp_HESS + ( dim_avg_n_Wrap(f1->SoilTemp(:,soil,{A:B},{C:D}),0)-273.15)\
                                 *Soil_thickness(soil)/sum(Soil_thickness)
      SoilMoist_HESS= SoilMoist_HESS + ( dim_avg_n_Wrap(f1->SoilMoist(:,soil,{A:B},{C:D}),0))\
                                 *Soil_thickness(soil)/sum(Soil_thickness)
    end do ; soil
    delete(f1)
end do ; year
Rain_LIS     = Rain_LIS/year_sum
Evap_LIS     = Evap_LIS/year_sum
TVeg_LIS     = TVeg_LIS/year_sum
ESoil_LIS    = ESoil_LIS/year_sum
ECanop_LIS   = ECanop_LIS/year_sum
Runoff_LIS   = Runoff_LIS/year_sum
Rnet_LIS     = Rnet_LIS/year_sum/12
Qle_LIS      = Qle_LIS/year_sum/12
Qh_LIS       = Qh_LIS/year_sum/12
Qg_LIS       = Qg_LIS/year_sum/12
SMsurf_LIS   = SMsurf_LIS/year_sum/12
LAI_LIS      = LAI_LIS/year_sum/12
SoilTemp_LIS = SoilTemp_LIS/year_sum/12
SoilMoist_LIS= SoilMoist_LIS/year_sum/12

Rain_HESS     = Rain_HESS/year_sum
Evap_HESS     = Evap_HESS/year_sum
TVeg_HESS     = TVeg_HESS/year_sum
ESoil_HESS    = ESoil_HESS/year_sum
ECanop_HESS   = ECanop_HESS/year_sum
Runoff_HESS   = Runoff_HESS/year_sum
Rnet_HESS     = Rnet_HESS/year_sum
Qle_HESS      = Qle_HESS/year_sum
Qh_HESS       = Qh_HESS/year_sum
Qg_HESS       = Qg_HESS/year_sum
LAI_HESS      = LAI_HESS/year_sum
SMsurf_HESS   = SMsurf_HESS/year_sum
SoilTemp_HESS = SoilTemp_HESS/year_sum
SoilMoist_HESS= SoilMoist_HESS/year_sum

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

; ================================ READ CLASS ==================================
file_CLASS = path_CLASS + "/CLASS_v1.1_"+year_s+".nc"
CL = addfile (file_CLASS,"r")

Rnet_CLASS= dim_avg_n_Wrap(CL->rs(:,{A:B},{C:D}),0)
Qle_CLASS = dim_avg_n_Wrap(CL->hfls(:,{A:B},{C:D}),0)
Qh_CLASS  = dim_avg_n_Wrap(CL->hfss(:,{A:B},{C:D}),0)
Qg_CLASS  = dim_avg_n_Wrap(CL->hfds(:,{A:B},{C:D}),0)
delete(CL)
do year = year_s+1,year_e

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



; ============================= plot water balance =============================
pic1 = "./plots/LIS_vs_obs_water_balance"
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed") ;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

; Rain
plots = new(2,graphic)
pres@txString = "Rain (mm/y)"
levels = fspan(0,1000,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,Rain_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,Rain_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)


; LAI
pres@txString = "LAI (mm/y)"
levels = fspan(0,6,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,LAI_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,LAI_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)
delete(plots)

; Evap
plots = new(3,graphic)
pres@txString = "Evap (mm/y)"
levels = fspan(0,1000,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,Evap_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,Evap_HESS,res)
res@tiMainString = "GLEAM"
plots(2) = gsn_csm_contour(wks1,Evap_GLEAM,res)
gsn_panel(wks1,(/plots/),(/1,3/),pres)

; TVeg
pres@txString = "TVeg (mm/y)"
levels = fspan(0,1000,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,TVeg_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,TVeg_HESS,res)
res@tiMainString = "GLEAM"
plots(2) = gsn_csm_contour(wks1,TVeg_GLEAM,res)
gsn_panel(wks1,(/plots/),(/1,3/),pres)

; ESoil
pres@txString = "ESoil (mm/y)"
levels = fspan(0,1000,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,ESoil_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,ESoil_HESS,res)
res@tiMainString = "GLEAM"
plots(2) = gsn_csm_contour(wks1,ESoil_GLEAM,res)
gsn_panel(wks1,(/plots/),(/1,3/),pres)

; ECanop
pres@txString = "ECanop (mm/y)"
levels = fspan(0,1000,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,ECanop_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,ECanop_HESS,res)
res@tiMainString = "GLEAM"
plots(2) = gsn_csm_contour(wks1,ECanop_GLEAM,res)
gsn_panel(wks1,(/plots/),(/1,3/),pres)

; Runoff
pres@txString = "Runoff (mm/y)"
levels = fspan(0,1000,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,Runoff_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,Runoff_HESS,res)
res@tiMainString = "LORA"
plots(2) = gsn_csm_contour(wks1,Runoff_LORA,res)
gsn_panel(wks1,(/plots/),(/1,3/),pres)

; SMsurf
pres@txString = "SMsurf (m3/m3)"
levels = fspan(0,0.5,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,SMsurf_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,SMsurf_HESS,res)
res@tiMainString = "GLEAM"
plots(2) = gsn_csm_contour(wks1,SMsurf_GLEAM,res)
gsn_panel(wks1,(/plots/),(/1,3/),pres)
delete(plots)

; ET/P
plots = new(2,graphic)
pres@txString = "Evap/Rain (-)"
levels = fspan(0,1.,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,Evap_LIS/Rain_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,Evap_HESS/Rain_HESS,res)
; res@tiMainString = "GLEAM"
; plots(2) = gsn_csm_contour(wks1,Evap_GLEAM/Rain_HESS,res)
; gsn_panel(wks1,(/plots/),(/1,3/),pres)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

; R/ET
pres@txString = "Runoff/Evap (-)"
levels = fspan(0,1.,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,Runoff_LIS/Evap_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,Runoff_HESS/Evap_HESS,res)
; res@tiMainString = "LORA/GLEAM"
; plots(2) = gsn_csm_contour(wks1,Runoff_LORA/Evap_GLEAM,res)
; gsn_panel(wks1,(/plots/),(/1,3/),pres)
gsn_panel(wks1,(/plots/),(/1,2/),pres)
delete(plots)
delete(levels)
delete(res@cnLevels)

; ============================ plot energy balance =============================
pic2 = "./plots/LIS_vs_obs_energy_balance"
wks2 = gsn_open_wks("pdf",pic2)
gsn_define_colormap(wks2,"WhiteBlueGreenYellowRed") ;"BlueWhiteOrangeRed") ;"BlueYellowRed")

; Rnet
plots = new(4,graphic)
pres@txString = "Rnet (W/m2)"
levels = fspan(-50,350,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks2,Rnet_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks2,Rnet_HESS,res)
res@tiMainString = "CLASS"
plots(2) = gsn_csm_contour(wks2,Rnet_CLASS,res)
res@tiMainString = "CERES"
plots(3) = gsn_csm_contour(wks2,Rnet_CERES,res)
gsn_panel(wks2,(/plots/),(/2,2/),pres)

; Qle
pres@txString = "Qle (W/m2)"
levels = fspan(-50,200,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks2,Qle_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks2,Qle_HESS,res)
res@tiMainString = "CLASS"
plots(2) = gsn_csm_contour(wks2,Qle_CLASS,res)
res@tiMainString = "DOLCE"
plots(3) = gsn_csm_contour(wks2,Qle_DOLCE,res)
gsn_panel(wks2,(/plots/),(/2,2/),pres)
delete(plots)

; Qh
plots = new(3,graphic)
pres@txString = "Qh (W/m2)"
levels = fspan(-50,200,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks2,Qh_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks2,Qh_HESS,res)
res@tiMainString = "CLASS"
plots(2) = gsn_csm_contour(wks2,Qh_CLASS,res)
gsn_panel(wks2,(/plots/),(/1,3/),pres)

; Qg
pres@txString = "Qg (W/m2)"
levels = fspan(-50,200,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks2,Qg_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks2,Qg_HESS,res)
res@tiMainString = "CLASS"
plots(2) = gsn_csm_contour(wks2,Qg_CLASS,res)
gsn_panel(wks2,(/plots/),(/1,3/),pres)


end
