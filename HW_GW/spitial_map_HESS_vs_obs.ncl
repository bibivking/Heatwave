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
case_name = (/"CTL","FREE_DRAIN", "PUMPx10","PUMPx100"/)
path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
case_sum  = dimsizes(case_name)
path      = new(case_sum,"string")
pic_message = "no-spinup"

do case_num = 0, case_sum -1
   path(case_num) = path_mum+case_name(case_num)+"/outputs-without-spinup"
end do

path_GLEAM= "/g/data/w35/Shared_data/Observations/Global_ET_products/GLEAM_v3_3/3_3a/yearly";
path_CLASS= "/g/data/w35/mm3972/data/CLASS_v1"
path_DOLCE= "/g/data/w35/mm3972/data/DOLCE" ; DOLCE_v1.0_2000.nc
path_LORA = "/g/data/w35/mm3972/data/LORA"  ; LORA_v1.0_2000.nc
path_CERES= "/g/data/w35/mm3972/data/CERES" ; CERES_EBAF-Surface_Ed2.8_Subset_200003-201502.nc

; ================================= Predefine ==================================
year_s = 2000
year_e = 2019
year_sum = year_e-year_s+1
row_sum = 2
col_sum = 3

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Soil_thickness1 = Soil_thickness
;(/0.005, 0.075, 0.154, 0.409, 1.085, 2.872/) ; LIS-CABLE


; ================================== OPTIONS ===================================
plot_type = "all"
scale     = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"

if scale .eq. "Global" then
  A = -90.
  B = 90.
  C = 0.
  D = 360.
else if scale .eq. "SE-AU" then
  A = -40.
  B = -28.
  C = 140.
  D = 154.
else if scale .eq. "AU" then
  A = -44.
  B = -10.
  C = 112.
  D = 154.
end if
end if
end if

filename  = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f         = addfile(filename,"r")

time_sum  = dimsizes(f->Rainf(:,0,0))
lat_sum   = dimsizes(f->Rainf(0,{A:B},0))
lon_sum   = dimsizes(f->Rainf(0,0,{C:D}))
lat_size  = dimsizes(f->y)
lon_size  = dimsizes(f->x)
CABLE_lat = f->y
CABLE_lon = f->x
default_val = -1.e+33
delete(f)

; Rain      = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
; Tair      = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
; Qair      = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
; Psurf     = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
; SWdown    = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
; LWdown    = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
; Wind      = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
; LAI       = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)

Evap      = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
TVeg      = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
ESoil     = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
ECanop    = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
Runoff    = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
Qrec      = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
Rnet      = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)

GWwb      = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
Qle       = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
Qh        = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
Qg        = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
SMsurf    = new((/case_sum, time_sum, lat_sum, lon_sum/),float, default_val)
SoilTemp  = new((/case_sum, time_sum, 6, lat_sum, lon_sum/),float, default_val)
SoilMoist = new((/case_sum, time_sum, 6, lat_sum, lon_sum/),float, default_val)

; Rain     = 0.
; Tair     = 0.
; Qair     = 0.
; Psurf    = 0.
; SWdown   = 0.
; LWdown   = 0.
; Wind     = 0.
; LAI      = 0.

Evap     = 0.
TVeg     = 0.
ESoil    = 0.
ECanop   = 0.
Runoff   = 0.
Qrec     = 0.
Rnet     = 0.
Qle      = 0.
Qh       = 0.
Qg       = 0.
GWwb     = 0.
SMsurf   = 0.
SoilTemp = 0.
SoilMoist= 0.

do case_num = 0, case_sum -1
   do year = year_s, year_e

      print("year = " +year)

      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if

      filename = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
      f        = addfile (filename,"r")

      ; Rain(case_num,year-year_s,:,:)     = dim_sum_n_Wrap(f->Rainf(:,{A:B},{C:D}),0)
      ; Tair(case_num,year-year_s,:,:)     = dim_avg_n_Wrap(f->Tair(:,{A:B},{C:D}),0)
      ; Qair(case_num,year-year_s,:,:)     = dim_avg_n_Wrap(f->Qair(:,{A:B},{C:D}),0)
      ; Psurf(case_num,year-year_s,:,:)    = dim_avg_n_Wrap(f->PSurf(:,{A:B},{C:D}),0)
      ; SWdown(case_num,year-year_s,:,:)   = dim_avg_n_Wrap(f->SWdown(:,{A:B},{C:D}),0)
      ; LWdown(case_num,year-year_s,:,:)   = dim_avg_n_Wrap(f->LWdown(:,{A:B},{C:D}),0)
      ; Wind(case_num,year-year_s,:,:)     = dim_avg_n_Wrap(f->Wind(:,{A:B},{C:D}),0)
      ; LAI(case_num,year-year_s,:,:)      = dim_avg_n_Wrap(f->LAI(:,{A:B},{C:D}),0)

      Evap(case_num,year-year_s,:,:)     = dim_sum_n_Wrap(f->Evap(:,{A:B},{C:D}),0)
      TVeg(case_num,year-year_s,:,:)     = dim_sum_n_Wrap(f->TVeg(:,{A:B},{C:D}),0)
      ESoil(case_num,year-year_s,:,:)    = dim_sum_n_Wrap(f->ESoil(:,{A:B},{C:D}),0)
      ECanop(case_num,year-year_s,:,:)   = dim_sum_n_Wrap(f->ECanop(:,{A:B},{C:D}),0)
      Runoff(case_num,year-year_s,:,:)   = dim_sum_n_Wrap(f->Qs(:,{A:B},{C:D}),0) \
                                           + dim_sum_n_Wrap(f->Qsb(:,{A:B},{C:D}),0)
      Qrec(case_num,year-year_s,:,:)     = dim_sum_n_Wrap(f->Qrecharge(:,{A:B},{C:D}),0)

      Rnet(case_num,year-year_s,:,:)     = dim_avg_n_Wrap(f->Rnet(:,{A:B},{C:D}),0)
      Qle(case_num,year-year_s,:,:)      = dim_avg_n_Wrap(f->Qle(:,{A:B},{C:D}),0)
      Qh(case_num,year-year_s,:,:)       = dim_avg_n_Wrap(f->Qh(:,{A:B},{C:D}),0)
      Qg(case_num,year-year_s,:,:)       = dim_avg_n_Wrap(f->Qg(:,{A:B},{C:D}),0)
      SMsurf(case_num,year-year_s,:,:)   = (dim_avg_n_Wrap(f->SoilMoist(:,0,{A:B},{C:D}),0)*0.022 \
                                            + dim_avg_n_Wrap(f->SoilMoist(:,1,{A:B},{C:D}),0)*0.058 \
                                            + dim_avg_n_Wrap(f->SoilMoist(:,2,{A:B},{C:D}),0)*(0.1 - 0.022 - 0.058)) \
                                            / 0.1
      GWwb(case_num,year-year_s,:,:)     = dim_avg_n_Wrap(f->GWMoist(:,{A:B},{C:D}),0)

      SoilTemp(case_num,year-year_s,:,:,:)  = dim_avg_n_Wrap(f->SoilTemp(:,:,{A:B},{C:D}),0)
      SoilMoist(case_num,year-year_s,:,:,:) = dim_avg_n_Wrap(f->SoilMoist(:,:,{A:B},{C:D}),0)

      delete(f)
   end do ; year
end do ; case_num

; Rain     = Rain*3600.*24.
; Tair     = Tair - 273.15

Evap     = Evap*3600.*24.
TVeg     = TVeg*3600.*24.
ESoil    = ESoil*3600.*24.
ECanop   = ECanop*3600.*24.
Runoff   = Runoff*3600.*24.
Qrec     = Qrec*3600.*24.
SoilTemp = SoilTemp - 273.15

; ================================ READ GLEAM ==================================
; use 2000-2019

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

Evap_GLEAM     = dim_avg_n_Wrap(E_tmp(time|20:38,lat|:, lon|:),0)
TVeg_GLEAM     = dim_avg_n_Wrap(T_tmp(time|20:38,lat|:, lon|:),0)
ESoil_GLEAM    = dim_avg_n_Wrap(Es_tmp(time|20:38,lat|:, lon|:),0)
ECanop_GLEAM   = dim_avg_n_Wrap(Ec_tmp(time|20:38,lat|:, lon|:),0)
SMsurf_GLEAM   = dim_avg_n_Wrap(SM_tmp(time|20:38,lat|:, lon|:),0)

;printVarSummary(Evap_GLEAM)
delete(G_ET)
delete(G_T)
delete(G_Es)
delete(G_Ec)
delete(G_SM)

; ================================ READ CLASS ==================================
; use 2003-2009

file_CLASS = path_CLASS + "/CLASS_v1.1_2003.nc"
CL         = addfile (file_CLASS,"r")
CLASS_lat  = doubletofloat(CL->lat({A:B}))
CLASS_lon  = doubletofloat(CL->lon({C:D}))

Rnet_CLASS= new((/dimsizes(CLASS_lat), dimsizes(CLASS_lon)/),float)
Qle_CLASS = new((/dimsizes(CLASS_lat), dimsizes(CLASS_lon)/),float)
Qh_CLASS  = new((/dimsizes(CLASS_lat), dimsizes(CLASS_lon)/),float)
Qg_CLASS  = new((/dimsizes(CLASS_lat), dimsizes(CLASS_lon)/),float)

Rnet_CLASS= 0.
Qle_CLASS = 0.
Qh_CLASS  = 0.
Qg_CLASS  = 0.

do year = 2003,2009
    file_CLASS = path_CLASS + "/CLASS_v1.1_"+year+".nc"
    CL = addfile (file_CLASS,"r")

    Rnet_CLASS = Rnet_CLASS + dim_avg_n_Wrap(CL->rs(:,{A:B},{C:D}),0)
    Qle_CLASS  = Qle_CLASS  + dim_avg_n_Wrap(CL->hfls(:,{A:B},{C:D}),0)
    Qh_CLASS   = Qh_CLASS   + dim_avg_n_Wrap(CL->hfss(:,{A:B},{C:D}),0)
    Qg_CLASS   = Qg_CLASS   + dim_avg_n_Wrap(CL->hfds(:,{A:B},{C:D}),0)
    delete(CL)
end do ;year

Rnet_CLASS = Rnet_CLASS/(2009-2003+1)
Qle_CLASS  = Qle_CLASS/(2009-2003+1)
Qh_CLASS   = Qh_CLASS/(2009-2003+1)
Qg_CLASS   = Qg_CLASS/(2009-2003+1)

Rnet_CLASS!0 = "lat"
Rnet_CLASS!1 = "lon"
Rnet_CLASS&lat = CLASS_lat
Rnet_CLASS&lon = CLASS_lon

Qle_CLASS!0 = "lat"
Qle_CLASS!1 = "lon"
Qle_CLASS&lat = CLASS_lat
Qle_CLASS&lon = CLASS_lon

Qh_CLASS!0 = "lat"
Qh_CLASS!1 = "lon"
Qh_CLASS&lat = CLASS_lat
Qh_CLASS&lon = CLASS_lon

Qg_CLASS!0 = "lat"
Qg_CLASS!1 = "lon"
Qg_CLASS&lat = CLASS_lat
Qg_CLASS&lon = CLASS_lon

delete(CLASS_lat)
delete(CLASS_lon)

; ================================ READ DOLCE ==================================
; ; 2000-2009
; 2003-2009

file_DOLCE = path_DOLCE +"/DOLCE_v1.0_"+year_s+".nc"
DOLCE      = addfile (file_DOLCE,"r")

DOLCE_lat  = doubletofloat(DOLCE->lat({A:B}))
DOLCE_lon  = doubletofloat(DOLCE->lon({C:D}))

Qle_DOLCE  = new((/dimsizes(DOLCE_lat), dimsizes(DOLCE_lon)/),float)
Qle_DOLCE  = 0.

do year = 2003,2009 ;2000, 2009
  file_DOLCE = path_DOLCE +"/DOLCE_v1.0_"+year+".nc"
  DOLCE      = addfile (file_DOLCE,"r")
  Qle_DOLCE  = Qle_DOLCE + dim_avg_n_Wrap(DOLCE->ET(:,{A:B},{C:D}),0)
  delete(DOLCE)
end do ;year
; Qle_DOLCE   = Qle_DOLCE/(2009-2000+1)
Qle_DOLCE   = Qle_DOLCE/(2009-2003+1)

Qle_DOLCE!0 = "lat"
Qle_DOLCE!1 = "lon"
Qle_DOLCE&lat = DOLCE_lat
Qle_DOLCE&lon = DOLCE_lon

delete(DOLCE_lat)
delete(DOLCE_lon)

; ================================ READ LORA ===================================
; 2000-2009

file_LORA   = path_LORA +"/LORA_v1.0_"+year_s+".nc"
L           = addfile (file_LORA,"r")

LORA_lat    = doubletofloat(L->lat({A:B}))
LORA_lon    = doubletofloat(L->lon({C:D}))

Runoff_LORA = new((/dimsizes(LORA_lat), dimsizes(LORA_lon)/),float)
Runoff_LORA = 0.

do year = 2000,2009
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

Runoff_LORA   = Runoff_LORA/(2009-2000+1)

Runoff_LORA!0 = "lat"
Runoff_LORA!1 = "lon"
Runoff_LORA&lat = LORA_lat
Runoff_LORA&lon = LORA_lon

delete(LORA_lat)
delete(LORA_lon)

; ================================ READ CERES ==================================
; 2001 -2018
file_CERES = path_CERES +"/CERES_EBAF-Surface_Ed4.1_Subset_200003-201911.nc"
CERES      = addfile (file_CERES,"r")
; 2001 -2018
; Rnet_CERES = dim_avg_n_Wrap(CERES->sfc_net_tot_all_mon(9:225,{A:B},{C:D}),0)
; 2003 -2009
Rnet_CERES = dim_avg_n_Wrap(CERES->sfc_net_tot_all_mon(9+24:225-9*12,{A:B},{C:D}),0)
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

; =============================== plot met input ===============================

;
; if (plot_type .eq. "met_force") .or. (plot_type .eq. "all") then
;
;     plots = new(dimsizes(case_name)+1,graphic)
;
;     pic0 = "./plots/HESS_vs_obs_met_"+pic_message+"_"+year_s+"-"+year_e
;     wks0 = gsn_open_wks("pdf",pic0)
;     gsn_define_colormap(wks0,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")
;
;     ; LAI
;     pres@txString = "LAI (-)"
;     res@cnLevels  = fspan(0,6,21)
;     do case_num = 0, case_sum-1
;       res@tiMainString = case_name(case_num)
;       plots(case_num) = gsn_csm_contour(wks1,dim_avg_n_Wrap(Rain(case_num,0:18,:,:),0),res)
;     end do
;     gsn_panel(wks1,(/plots/),(/row_sum,col_sum/),pres)
;     delete(res@cnLevels)
;
;     ; Rain
;     pres@txString = "Rain (mm/y)"
;     res@cnLevels  = fspan(0,1000,21)
;     do case_num = 0, case_sum-1
;       res@tiMainString = case_name(case_num)
;       plots(case_num) = gsn_csm_contour(wks0,dim_avg_n_Wrap(Rain(case_num,:,:,:),0),res)
;     end do
;     gsn_panel(wks0,(/plots/),(/row_sum,col_sum/),pres)
;     delete(res@cnLevels)
;
;     ; Tair
;     pres@txString = "Tair (C)"
;     res@cnLevels  = fspan(0,30,16)
;     res@cnFillPalette = "WhiteBlueGreenYellowRed"
;     res@tiMainString = case_name(0)
;     plots(0) = gsn_csm_contour(wks0,dim_avg_n_Wrap(Tair(0,:,:,:),0),res)
;     delete(res@cnLevels)
;     do case_num = 1, case_sum-1
;       res@cnFillPalette = "ViBlGrWhYeOrRe"
;       res@tiMainString = case_name(case_num) +" - "+ case_name(0)
;       res@cnLevels     = fspan(-5.,5.,10)
;       plots(case_num)  = gsn_csm_contour(wks0,dim_avg_n_Wrap(Tair(case_num,:,:,:),0)-dim_avg_n_Wrap(Tair(0,:,:,:),0),res)
;     end do
;     gsn_panel(wks0,(/plots/),(/row_sum,col_sum/),pres)
;     delete(res@cnLevels)
;
;     ; Qair
;     pres@txString =  "Qair (%)"
;     res@cnLevels  = fspan(0,1,21)
;     res@cnFillPalette = "WhiteBlueGreenYellowRed"
;     res@tiMainString = case_name(0)
;     plots(0) = gsn_csm_contour(wks0,dim_avg_n_Wrap(Qair(0,:,:,:),0)*100.,res)
;     delete(res@cnLevels)
;     do case_num = 1, case_sum-1
;       res@cnFillPalette = "ViBlGrWhYeOrRe"
;       res@tiMainString = case_name(case_num) +" - "+ case_name(0)
;       res@cnLevels     = fspan(-0.5,0.5,21)
;       plots(case_num)  = gsn_csm_contour(wks0, dim_avg_n_Wrap(Qair(case_num,:,:,:),0)*100. - \
;                                                dim_avg_n_Wrap(Qair(0,:,:,:),0)*100.,res)
;     end do
;     gsn_panel(wks0,(/plots/),(/row_sum,col_sum/),pres)
;     delete(res@cnLevels)
;
;
;     ; Psurf
;     pres@txString = "Psurf (hPa)"
;     res@cnLevels  = fspan(900,1200,31)
;     res@cnFillPalette = "WhiteBlueGreenYellowRed"
;     res@tiMainString = case_name(0)
;     plots(0) = gsn_csm_contour(wks0,dim_avg_n_Wrap(Psurf(0,:,:,:),0)*100.,res)
;     delete(res@cnLevels)
;     do case_num = 1, case_sum-1
;       res@cnFillPalette = "ViBlGrWhYeOrRe"
;       res@tiMainString = case_name(case_num) +" - "+ case_name(0)
;       res@cnLevels     = fspan(-100,100,21)
;       plots(case_num)  = gsn_csm_contour(wks0, dim_avg_n_Wrap(Psurf(case_num,:,:,:),0)*100. - \
;                                                dim_avg_n_Wrap(Psurf(0,:,:,:),0)*100.,res)
;     end do
;     gsn_panel(wks0,(/plots/),(/row_sum,col_sum/),pres)
;     delete(res@cnLevels)
;
;     ; SWdown
;     pres@txString = "SWdown (W/m^2)"
;     res@cnLevels  = fspan(0,500,51)
;     res@cnFillPalette = "WhiteBlueGreenYellowRed"
;     res@tiMainString = case_name(0)
;     plots(0) = gsn_csm_contour(wks0,dim_avg_n_Wrap(SWdown(0,:,:,:),0),res)
;     delete(res@cnLevels)
;     do case_num = 1, case_sum-1
;       res@cnFillPalette = "ViBlGrWhYeOrRe"
;       res@tiMainString = case_name(case_num) +" - "+ case_name(0)
;       res@cnLevels     = fspan(-100,100,41)
;       plots(case_num)  = gsn_csm_contour(wks0, dim_avg_n_Wrap(SWdown(case_num,:,:,:),0) - \
;                                                dim_avg_n_Wrap(SWdown(0,:,:,:),0),res)
;     end do
;     gsn_panel(wks0,(/plots/),(/row_sum,col_sum/),pres)
;     delete(res@cnLevels)
;
;     ; LWdown
;     pres@txString = "LWdown (W/m^2)"
;     res@cnLevels  = fspan(0,500,51)
;     res@cnFillPalette = "WhiteBlueGreenYellowRed"
;     res@tiMainString = case_name(0)
;     plots(0) = gsn_csm_contour(wks0,dim_avg_n_Wrap(LWdown(0,:,:,:),0),res)
;     delete(res@cnLevels)
;     do case_num = 1, case_sum-1
;       res@cnFillPalette = "ViBlGrWhYeOrRe"
;       res@tiMainString = case_name(case_num) +" - "+ case_name(0)
;       res@cnLevels     = fspan(-100,100,41)
;       plots(case_num)  = gsn_csm_contour(wks0, dim_avg_n_Wrap(LWdown(case_num,:,:,:),0) - \
;                                                dim_avg_n_Wrap(LWdown(0,:,:,:),0),res)
;     end do
;     gsn_panel(wks0,(/plots/),(/row_sum,col_sum/),pres)
;     delete(res@cnLevels)
;
;     ; Wind
;     pres@txString = "Wind (m/s)"
;     res@cnLevels  = fspan(0,10,21)
;     res@cnFillPalette = "WhiteBlueGreenYellowRed"
;     res@tiMainString = case_name(0)
;     plots(0) = gsn_csm_contour(wks0,dim_avg_n_Wrap(Wind(0,:,:,:),0),res)
;     delete(res@cnLevels)
;     do case_num = 1, case_sum-1
;       res@cnFillPalette = "ViBlGrWhYeOrRe"
;       res@tiMainString = case_name(case_num) +" - "+ case_name(0)
;       res@cnLevels     = fspan(-3,3,21)
;       plots(case_num)  = gsn_csm_contour(wks0, dim_avg_n_Wrap(Wind(case_num,:,:,:),0) - \
;                                                dim_avg_n_Wrap(Wind(0,:,:,:),0),res)
;     end do
;     gsn_panel(wks0,(/plots/),(/row_sum,col_sum/),pres)
;     delete(res@cnLevels)
;     delete(plots)
; end if

; ============================== plot water balance ============================
if (plot_type .eq. "water_bal") .or. (plot_type .eq. "all") then

    plots = new(dimsizes(case_name)+1,graphic)

    pic1 = "./plots/HESS_vs_obs_water_balance_"+pic_message+"_"+year_s+"-"+year_e
    wks1 = gsn_open_wks("pdf",pic1)
    gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    ; Evap
    pres@txString = "Evap (mm/y)"
    res@cnLevels  = fspan(0,1000,21)
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks1, dim_avg_n_Wrap(Evap(case_num,0:18,:,:),0),res)
    end do
    res@tiMainString = "GLEAM"
    plots(case_sum) = gsn_csm_contour(wks1,Evap_GLEAM,res)
    gsn_panel(wks1,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)

    ; TVeg
    pres@txString = "TVeg (mm/y)"
    res@cnLevels  = fspan(0,1000,21)
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks1, dim_avg_n_Wrap(TVeg(case_num,0:18,:,:),0),res)
    end do
    res@tiMainString = "GLEAM"
    plots(case_sum)  = gsn_csm_contour(wks1,TVeg_GLEAM,res)
    gsn_panel(wks1,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)

    ; ESoil
    pres@txString = "ESoil (mm/y)"
    res@cnLevels  = fspan(0,1000,21)
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks1, dim_avg_n_Wrap(ESoil(case_num,0:18,:,:),0),res)
    end do
    res@tiMainString = "GLEAM"
    plots(case_sum)  = gsn_csm_contour(wks1,ESoil_GLEAM,res)
    gsn_panel(wks1,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)

    ; ECanop
    pres@txString = "ECanop (mm/y)"
    res@cnLevels  = fspan(0,1000,21)
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks1, dim_avg_n_Wrap(ECanop(case_num,0:18,:,:),0),res)
    end do
    res@tiMainString = "GLEAM"
    plots(case_sum)  = gsn_csm_contour(wks1,ECanop_GLEAM,res)
    gsn_panel(wks1,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)

    ; Runoff
    pres@txString = "Runoff (mm/y)"
    res@cnLevels  = fspan(0,1000,21)
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks1, dim_avg_n_Wrap(Runoff(case_num,0:9,:,:),0),res)
    end do
    res@tiMainString = "LORA"
    plots(case_sum)  = gsn_csm_contour(wks1,Runoff_LORA,res)
    gsn_panel(wks1,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)
    delete(plots(case_sum))

    ; Recharge
    pres@txString = "Recharge (mm/y)"
    res@cnLevels  = fspan(-200,200,21)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks1, dim_avg_n_Wrap(Qrec(case_num,0:18,:,:),0),res)
    end do
    gsn_panel(wks1,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)

    ; SMsurf
    pres@txString = "SMsurf (m3/m3)"
    res@cnLevels  = fspan(0,0.5,21)
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks1, dim_avg_n_Wrap(SMsurf(case_num,0:18,:,:),0),res)
    end do
    res@tiMainString = "GLEAM"
    plots(case_sum)  = gsn_csm_contour(wks1,SMsurf_GLEAM,res)
    gsn_panel(wks1,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)
    delete(plots(case_sum))

    ; Soil Moist
    do soil = 0,5
        pres@txString = "Soil Moisture in "+(soil+1)+" layer (m3/m3)"
        res@cnLevels  = fspan(0,0.5,21)
        res@cnFillPalette = "WhiteBlueGreenYellowRed"
        do case_num = 0, case_sum-1
          res@tiMainString = case_name(case_num)
          plots(case_num)  = gsn_csm_contour(wks1, dim_avg_n_Wrap(SoilMoist(case_num,0:18,soil,:,:),0),res)
        end do
        gsn_panel(wks1,(/plots/),(/row_sum,col_sum/),pres)
        delete(res@cnLevels)
    end do

    ; GWwb
    pres@txString = "Aquifer Moisture (m3/m3)"
    res@cnLevels  = fspan(0,0.5,21)
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks1, dim_avg_n_Wrap(GWwb(case_num,0:18,:,:),0),res)
    end do
    gsn_panel(wks1,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)

    ; Soil Temp
    do soil = 0,5
        pres@txString = "Soil Temperature in "+(soil+1)+" layer (C)"
        res@cnLevels  = fspan(0,40,21)
        res@cnFillPalette = "WhiteBlueGreenYellowRed"
        do case_num = 0, case_sum-1
          res@tiMainString = case_name(case_num)
          plots(case_num)  = gsn_csm_contour(wks1, dim_avg_n_Wrap(SoilTemp(case_num,0:18,soil,:,:),0),res)
        end do
        gsn_panel(wks1,(/plots/),(/row_sum,col_sum/),pres)
        delete(res@cnLevels)
    end do

    ; ; ET/P
    ; plots = new(2,graphic)
    ; pres@txString = "Evap/Rain (-)"
    ; levels = fspan(0,1.,21)
    ; res@cnLevels = levels
    ; res@cnFillPalette = "WhiteBlueGreenYellowRed"
    ; res@tiMainString = "CTL"
    ; plots(0) = gsn_csm_contour(wks1,Evap/Rain,res)
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
    ; plots(0) = gsn_csm_contour(wks1,Runoff/Evap,res)
    ; res@tiMainString = "SEN"
    ; plots(1) = gsn_csm_contour(wks1,Runoff_SEN/Evap_SEN,res)
    ; ; res@tiMainString = "LORA/GLEAM"
    ; ; plots(2) = gsn_csm_contour(wks1,Runoff_LORA/Evap_GLEAM,res)
    ; ; gsn_panel(wks1,(/plots/),(/1,3/),pres)
    ; gsn_panel(wks1,(/plots/),(/1,2/),pres)
    ; delete(plots)
    ; delete(levels)
    ; delete(res@cnLevels)
    delete(plots)
end if

; ============================ plot energy balance =============================
if (plot_type .eq. "energy_bal") .or. (plot_type .eq. "all") then

    plots = new(dimsizes(case_name)+2,graphic)

    pic2 = "./plots/HESS_vs_obs_energy_balance_"+pic_message+"_"+year_s+"-"+year_e
    wks2 = gsn_open_wks("pdf",pic2)
    gsn_define_colormap(wks2,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    ; Rnet
    pres@txString = "Rnet (W/m2)"
    res@cnLevels  = fspan(-50,200,51)
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks2, dim_avg_n_Wrap(Rnet(case_num,3:9,:,:),0),res)
    end do
    res@tiMainString = "CERES"
    plots(case_sum)  = gsn_csm_contour(wks2,Rnet_CERES,res)
    res@tiMainString = "CLASS"
    plots(case_sum+1)  = gsn_csm_contour(wks2,Rnet_CLASS,res)
    gsn_panel(wks2,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)

    ; Qle
    pres@txString = "Qle (W/m2)"
    res@cnLevels  = fspan(-50,200,51)
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks2, dim_avg_n_Wrap(Qle(case_num,3:9,:,:),0),res)
    end do
    res@tiMainString = "CLASS"
    plots(case_sum)  = gsn_csm_contour(wks2,Qle_CLASS,res)
    res@tiMainString = "DOLCE"
    plots(case_sum+1)= gsn_csm_contour(wks2,Qle_DOLCE,res)
    gsn_panel(wks2,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)
    delete(plots(case_sum+1))

    ; Qh
    pres@txString = "Qh (W/m2)"
    res@cnLevels  = fspan(-50,200,51)
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks2, dim_avg_n_Wrap(Qh(case_num,3:9,:,:),0),res)
    end do
    res@tiMainString = "CLASS"
    plots(case_sum)  = gsn_csm_contour(wks2,Qh_CLASS,res)
    gsn_panel(wks2,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)

    ; Qg
    pres@txString = "Qg (W/m2)"
    res@cnLevels  = fspan(-50,200,51)
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0, case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num)  = gsn_csm_contour(wks2, dim_avg_n_Wrap(Qg(case_num,3:9,:,:),0),res)
    end do
    res@tiMainString = "CLASS"
    plots(case_sum)  = gsn_csm_contour(wks2,Qg_CLASS,res)
    gsn_panel(wks2,(/plots/),(/row_sum,col_sum/),pres)
    delete(res@cnLevels)
end if
exit
end
