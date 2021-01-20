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


; ================================= Predefine ==================================
pic_name = "CTL_vs_PUMP"
case_name = (/"CTL","PUMP","PUMPx10","PUMPx100"/);,"PUMP_IRRG","PUMP_IRRGx10","PUMP_IRRGx100"/)
case_sum  = dimsizes(case_name)

year_s = 2000
year_e = 2019
year_sum = year_e-year_s+1

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

; =================================== Path =====================================
path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"+case_name(case_num)+"/outputs"
end do

path_GLEAM= "/g/data/w35/Shared_data/Observations/Global_ET_products/GLEAM_v3_3/3_3a/yearly";
path_CLASS= "/g/data/w35/mm3972/data/CLASS_v1"
path_DOLCE= "/g/data/w35/mm3972/data/DOLCE" ; DOLCE_v1.0_2000.nc
path_LORA = "/g/data/w35/mm3972/data/LORA"  ; LORA_v1.0_2000.nc
path_CERES= "/g/data/w35/mm3972/data/CERES" ; CERES_EBAF-Surface_Ed2.8_Subset_200003-201502.nc

; ================================== OPTIONS ===================================
plot_type = "all"
scale = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"

filename = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
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

lat_sum = dimsizes(f->Rainf(0,{A:B},0))
lon_sum = dimsizes(f->Rainf(0,0,{C:D}))
delete(f)

filename1 = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f1        = addfile (filename1,"r")

Rain     = new((/case_sum,lat_sum,lon_sum/),float,f1->Rainf@_FillValue)
Tair     = new((/case_sum,lat_sum,lon_sum/),float,f1->Tair@_FillValue)
LAI      = new((/case_sum,lat_sum,lon_sum/),float,f1->LAI@_FillValue)

Evap     = new((/case_sum,lat_sum,lon_sum/),float,f1->Evap@_FillValue)
TVeg     = new((/case_sum,lat_sum,lon_sum/),float,f1->TVeg@_FillValue)
ESoil    = new((/case_sum,lat_sum,lon_sum/),float,f1->ESoil@_FillValue)
ECanop   = new((/case_sum,lat_sum,lon_sum/),float,f1->ECanop@_FillValue)
Runoff   = new((/case_sum,lat_sum,lon_sum/),float,f1->Qs@_FillValue)
Qrec     = new((/case_sum,lat_sum,lon_sum/),float,f1->Qrecharge@_FillValue)
WatTable = new((/case_sum,lat_sum,lon_sum/),float,f1->WatTable@_FillValue)
Rnet     = new((/case_sum,lat_sum,lon_sum/),float,f1->Rnet@_FillValue)
Qle      = new((/case_sum,lat_sum,lon_sum/),float,f1->Qle@_FillValue)
Qh       = new((/case_sum,lat_sum,lon_sum/),float,f1->Qh@_FillValue)
Qg       = new((/case_sum,lat_sum,lon_sum/),float,f1->Qg@_FillValue)
SMsurf   = new((/case_sum,lat_sum,lon_sum/),float,f1->SoilMoist@_FillValue)
SoilTemp = new((/case_sum,6,lat_sum,lon_sum/),float,f1->SoilTemp@_FillValue)
SoilMoist= new((/case_sum,6,lat_sum,lon_sum/),float,f1->SoilMoist@_FillValue)
GWwb     = new((/case_sum,lat_sum,lon_sum/),float,f1->GWMoist@_FillValue)

Rain     = 0.
Tair     = 0.
Evap     = 0.
TVeg     = 0.
ESoil    = 0.
ECanop   = 0.
Runoff   = 0.
Qrec     = 0.
WatTable = 0.
Rnet     = 0.
Qle      = 0.
Qh       = 0.
Qg       = 0.
SoilTemp = 0.
SoilMoist= 0.
SMsurf   = 0.
LAI      = 0.
GWwb     = 0.
delete(f1)

do case_num = 0,case_sum-1
  do year = year_s,year_e
      print("year = " +year)

      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if

      filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
      f1 = addfile (filename1,"r")

      Rain(case_num,:,:)     = Rain(case_num,:,:)   + dim_avg_n_Wrap(f1->Rainf(:,{A:B},{C:D}),0)*3600*24
      Evap(case_num,:,:)     = Evap(case_num,:,:)   + dim_avg_n_Wrap(f1->Evap(:,{A:B},{C:D}),0)*3600*24
      TVeg(case_num,:,:)     = TVeg(case_num,:,:)   + dim_avg_n_Wrap(f1->TVeg(:,{A:B},{C:D}),0)*3600*24
      ESoil(case_num,:,:)    = ESoil(case_num,:,:)  + dim_avg_n_Wrap(f1->ESoil(:,{A:B},{C:D}),0)*3600*24
      ECanop(case_num,:,:)   = ECanop(case_num,:,:) + dim_avg_n_Wrap(f1->ECanop(:,{A:B},{C:D}),0)*3600*24
      Runoff(case_num,:,:)   = Runoff(case_num,:,:) + dim_avg_n_Wrap(f1->Qs(:,{A:B},{C:D}),0)*3600*24 \
                                                    + dim_avg_n_Wrap(f1->Qsb(:,{A:B},{C:D}),0)*3600*24
      Qrec(case_num,:,:)     = Qrec(case_num,:,:) + dim_avg_n_Wrap(f1->Qrecharge(:,{A:B},{C:D}),0)*3600*24
      WatTable(case_num,:,:) = WatTable(case_num,:,:) + dim_avg_n_Wrap(f1->WatTable(:,{A:B},{C:D}),0)

      Tair(case_num,:,:)     = Tair(case_num,:,:)   + dim_avg_n_Wrap(f1->Tair(:,{A:B},{C:D}),0)-273.15
      LAI(case_num,:,:)      = LAI(case_num,:,:)    + dim_avg_n_Wrap(f1->LAI(:,{A:B},{C:D}),0)
      Rnet(case_num,:,:)     = Rnet(case_num,:,:)   + dim_avg_n_Wrap(f1->Rnet(:,{A:B},{C:D}),0)
      Qle(case_num,:,:)      = Qle(case_num,:,:)    + dim_avg_n_Wrap(f1->Qle(:,{A:B},{C:D}),0)
      Qh(case_num,:,:)       = Qh(case_num,:,:)     + dim_avg_n_Wrap(f1->Qh(:,{A:B},{C:D}),0)
      Qg(case_num,:,:)       = Qg(case_num,:,:)     + dim_avg_n_Wrap(f1->Qg(:,{A:B},{C:D}),0)
      SMsurf(case_num,:,:)   = SMsurf(case_num,:,:) + ( dim_avg_n_Wrap(f1->SoilMoist(:,0,{A:B},{C:D}),0)*0.022 \
                          + dim_avg_n_Wrap(f1->SoilMoist(:,1,{A:B},{C:D}),0)*0.058 \
                          + dim_avg_n_Wrap(f1->SoilMoist(:,2,{A:B},{C:D}),0)*(0.1 - 0.022 - 0.058)) \
                          / 0.1
      GWwb(case_num,:,:)     = GWwb(case_num,:,:)  + dim_avg_n_Wrap(f1->GWMoist(:,{A:B},{C:D}),0)

      SoilTemp(case_num,:,:,:) = SoilTemp(case_num,:,:,:)  + dim_avg_n_Wrap(f1->SoilTemp(:,:,{A:B},{C:D}),0)
      SoilMoist(case_num,:,:,:)= SoilMoist(case_num,:,:,:) + dim_avg_n_Wrap(f1->SoilMoist(:,:,{A:B},{C:D}),0)

      delete(f1)
    end do ; year
end do ; case_num
Rain     = Rain/year_sum
Evap     = Evap/year_sum
TVeg     = TVeg/year_sum
ESoil    = ESoil/year_sum
ECanop   = ECanop/year_sum
Runoff   = Runoff/year_sum
Qrec     = Qrec/year_sum
WatTable = WatTable/year_sum
Tair     = Tair/year_sum
Rnet     = Rnet/year_sum
Qle      = Qle/year_sum
Qh       = Qh/year_sum
Qg       = Qg/year_sum
LAI      = LAI/year_sum
SMsurf   = SMsurf/year_sum
SoilTemp = SoilTemp/year_sum
SoilMoist= SoilMoist/year_sum
GWwb     = GWwb/year_sum

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
file_CLASS = path_CLASS + "/CLASS_v1.1_"+year_s+"_SE_Aus.nc"
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

    file_CLASS = path_CLASS + "/CLASS_v1.1_"+year+"_SE_Aus.nc"
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
file_DOLCE = path_DOLCE +"/DOLCE_v1.0_"+year_s+"_SE_Aus.nc"
;print(file_DOLCE)
DOLCE      = addfile (file_DOLCE,"r")
Qle_DOLCE  = dim_avg_n_Wrap(DOLCE->ET(:,{A:B},{C:D}),0)

do year = year_s+1,year_e
  file_DOLCE = path_DOLCE +"/DOLCE_v1.0_"+year+"_SE_Aus.nc"
  DOLCE      = addfile (file_DOLCE,"r")
  Qle_DOLCE  = Qle_DOLCE + dim_avg_n_Wrap(DOLCE->ET(:,{A:B},{C:D}),0)
  delete(DOLCE)
end do ;year
Qle_DOLCE   = Qle_DOLCE/year_sum

; ================================ READ LORA ===================================
file_LORA   = path_LORA +"/LORA_v1.0_"+year_s+"_SE_Aus.nc"
L           = addfile (file_LORA,"r")
Runoff_LORA = dim_avg_n_Wrap(L->mrro(:,{A:B},{C:D}),0)
Runoff_LORA = 0.

do year = year_s,year_e
  if ( mod(year,4) .eq. 0) then
     dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
  else
     dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  end if

  file_LORA   = path_LORA +"/LORA_v1.0_"+year+"_SE_Aus.nc"
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

  ; ============================= plot water balance ===========================
  if (plot_type .eq. "water_bal") .or. (plot_type .eq. "all") then
    pic1 = "./plots/HESS_version_compare_water_balance_"+year_s+"-"+year_e+"_"+pic_name
    wks1 = gsn_open_wks("pdf",pic1)
    gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")
    plots = new(case_sum+1,graphic)

    ; LAI
    pres@txString = "LAI"
    levels = fspan(0,6,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks1,LAI(case_num,:,:),res)
    end do
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)

    ; Evap
    pres@txString = "Evap (mm/d)"
    levels = fspan(0,10,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks1,Evap(case_num,:,:),res)
    end do
    res@tiMainString = "GLEAM"
    plots(case_sum) = gsn_csm_contour(wks1,Evap_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)

    ; TVeg
    pres@txString = "TVeg (mm/d)"
    levels = fspan(0,10,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks1,TVeg(case_num,:,:),res)
    end do
    res@tiMainString = "GLEAM"
    plots(case_sum) = gsn_csm_contour(wks1,TVeg_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)


    ; ESoil
    pres@txString = "ESoil (mm/d)"
    levels = fspan(0,10,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks1,ESoil(case_num,:,:),res)
    end do
    res@tiMainString = "GLEAM"
    plots(case_sum) = gsn_csm_contour(wks1,ESoil_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)

    ; ECanop
    pres@txString = "ECanop (mm/d)"
    levels = fspan(0,10,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks1,ECanop(case_num,:,:),res)
    end do
    res@tiMainString = "GLEAM"
    plots(case_sum) = gsn_csm_contour(wks1,ECanop_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)

    ; Runoff
    pres@txString = "Runoff (mm/d)"
    levels = fspan(0,10,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks1,ECanop(case_num,:,:),res)
    end do
    res@tiMainString = "LORA"
    plots(case_sum) = gsn_csm_contour(wks1,Runoff_LORA,res)
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)

    ; Recharge
    pres@txString = "Recharge (mm/d)"
    levels = fspan(0,10,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks1,Qrec(case_num,:,:),res)
    end do
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)

    ; WatTable
    pres@txString = "WatTable (mm/d)"
    levels = fspan(0,100,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks1,WatTable(case_num,:,:),res)
    end do
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)

    ; SMsurf
    pres@txString = "SMsurf (mm/d)"
    levels = fspan(0,0.5,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks1,SMsurf(case_num,:,:),res)
    end do
    res@tiMainString = "GLEAM"
    plots(case_sum) = gsn_csm_contour(wks1,SMsurf_GLEAM,res)
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)

    ; Soil Moist
    do soil = 0,5
        pres@txString = "Soil Moisture of layer-"+soil+1+" (m3/m3)"
        levels = fspan(0,0.5,21)
        res@cnLevels = levels
        res@cnFillPalette = "WhiteBlueGreenYellowRed"
        do case_num = 0,case_sum-1
          res@tiMainString = case_name(case_num)
          plots(case_num) = gsn_csm_contour(wks1,SoilMoist(case_num,soil,:,:),res)
        end do
        delete(levels)
        delete(res@cnLevels)
        gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)
    end do

    pres@txString = "Aquifer Moisture (m3/m3)"
    levels = fspan(0,0.5,21)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks1,GWwb(case_num,:,:),res)
    end do
    delete(levels)
    delete(res@cnLevels)
    gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)

    ; Soil Temp
    do soil = 0,5
        pres@txString = "Soil Temperature of layer-"+soil+1+" (C)"
        levels = fspan(0,45,46)
        res@cnLevels = levels
        res@cnFillPalette = "WhiteBlueGreenYellowRed"
        do case_num = 0,case_sum-1
          res@tiMainString = case_name(case_num)
          plots(case_num) = gsn_csm_contour(wks1,SoilTemp(case_num,soil,:,:),res)
        end do
        delete(levels)
        delete(res@cnLevels)
        gsn_panel(wks1,(/plots/),(/ceil(case_sum/2),2/),pres)
    end do

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
    ; delete(levels)
    ; delete(res@cnLevels)

    delete(plots)
end if

; ============================ plot energy balance =============================
  if (plot_type .eq. "energy_bal") .or. (plot_type .eq. "all") then
    pic2 = "./plots/HESS_version_compare_energy_balance_"+year_s+"-"+year_e+"_"+pic_name
    wks2 = gsn_open_wks("pdf",pic2)
    gsn_define_colormap(wks2,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    plots = new(case_sum+1,graphic)

    ; Rnet
    pres@txString = "Rnet (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks2,Rnet(case_num,:,:),res)
    end do
    res@tiMainString = "CERES"
    plots(case_sum) = gsn_csm_contour(wks2,Rnet_CERES,res)
    gsn_panel(wks2,(/plots/),(/ceil(case_sum/2),2/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; Qle
    pres@txString = "Qle (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks2,Qle(case_num,:,:),res)
    end do
    res@tiMainString = "CLASS"
    plots(case_num) = gsn_csm_contour(wks2,Qle_CLASS,res)
    gsn_panel(wks2,(/plots/),(/ceil(case_sum/2),2/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; Qh
    pres@txString = "Qh (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks2,Qh(case_num,:,:),res)
    end do
    res@tiMainString = "CLASS"
    plots(case_num) = gsn_csm_contour(wks2,Qh_CLASS,res)
    gsn_panel(wks2,(/plots/),(/ceil(case_sum/2),2/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; Qg
    pres@txString = "Qg (W/m2)"
    levels = fspan(-50,200,51)
    res@cnLevels = levels
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    do case_num = 0,case_sum-1
      res@tiMainString = case_name(case_num)
      plots(case_num) = gsn_csm_contour(wks2,Qg(case_num,:,:),res)
    end do
    res@tiMainString = "CLASS"
    plots(case_num) = gsn_csm_contour(wks2,Qg_CLASS,res)
    gsn_panel(wks2,(/plots/),(/ceil(case_sum/2),2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end if

end
