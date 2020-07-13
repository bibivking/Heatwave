;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

path = "/g/data/w35/Shared_data/Observations/AWAP_all_variables/daily/"

ftmean = path + "tmean/AWAP_daily_tmean_1970_2019.nc"; tmean: minimum temperature (C) ???
ftmin  = path + "tmin/AWAP_daily_tmin_1970_2019.nc" ;tmin: minimum temperature (C)
ftmax  = path + "tmax/AWAP_daily_tmax_1970_2019.nc" ; tmax: maximum temperature (C)
frain  = path + "rain/AWAP_daily_rain_1970_2019.nc" ; rain: total precipitation (mm/d)
frad   = path + "rad/AWAP_daily_rad_1990_2019.nc" ; rad: incoming shortwave radiation (W/m2)
fwind  = path + "windspeed/AWAP_daily_windspeed_1970_2019.nc" ; windspeed: wind speed (m/s)
fvpd15 = path + "vph15/AWAP_daily_vph15_1970_2019.nc" ;
fvpd09 = path + "vph09/AWAP_daily_vph09_1970_2019.nc"


var_name= (/\
           "Evap","TVeg","ESoil",\ ;"kg/m^2/s"
           "Qs","Qsb",\      ;"kg/m^2/s"
           "Qle","Qh","Qg",\      ;"W/m^2"
           "Rnet","Albedo",\
           "RadT","VegT",\   ; "K"
           "SnowT","BaresoilT",\
           "GPP","NPP",\     ;"umol/m^2/s"
           "SoilTemp",\
           "SoilMoist",\     ;"m^3/m^3"  Average layer soil moisture
           "SoilMoistIce",\
           "GWMoist",\       ;"mm3/mm3" Aquifer mositure content, Groundwater_off doesn't have
           "WatTable"\       ;"m"
           /)
           ;"LAI" ; the same in Groundwater_on and Groundwater_off

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

;______________________________ Data _________________________________
year_s = 2003
year_e = 2009
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
case_sum = dimsizes(case_name)

scale = "Global" ;"Global" ; "50S-N" ; "SE-AU" ; "AU"

CABLE = new((/10,360,720/),float)
CABLE = 0.0
CLASS = new((/2,360,720/),float)
CLASS = 0.0
DOLCE = new((/360,720/),float)
DOLCE = 0.0

do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    filename = case_name(0)+"/cable_out_"+year+".nc"
    f = addfile (filename,"r")

    lat = f->Qle&y
    lon = f->Qle&x

    CABLE!0 = "var"
    CABLE!1 = "lat"
    CABLE!2 = "lon"
    CABLE&lat = lat
    CABLE&lon = lon

    CABLE@_FillValue   = f->Qle@_FillValue

    do var_num = 0,9
       if (var_num .le. 4) then
         do month = 0,11
            CABLE(var_num,:,:) = CABLE(var_num,:,:) + f->$var_name(var_num)$(month,:,:)*3600*24*dom(month)
         end do
       else
          CABLE(var_num,:,:) = CABLE(var_num,:,:) + dim_avg_n_Wrap(f->$var_name(var_num)$,0)
       end if
    end do ;var_num
    ; __________________ CLASS fluxes ______________________
    filename1 = case_name(1)+"/CLASS_v1.1_"+year+".nc"
    f1 = addfile (filename1,"r")

    CLASS(0,:,:) = CLASS(0,:,:) + dim_avg_n_Wrap(lonFlip(f1->hfls),0)*1.
    CLASS(1,:,:) = CLASS(1,:,:) + dim_avg_n_Wrap(lonFlip(f1->hfss),0)*1.

    ; ___________________ DOLCE Latent _________________________
    filename2 = case_name(2)+"/DOLCE_v1.0_"+year+".nc"
    f2 = addfile (filename2,"r")
    DOLCE = DOLCE + dim_avg_n_Wrap(lonFlip(f2->ET),0)*1.0

end do ;year

CLASS!0 = "var"
CLASS!1 = "lat"
CLASS!2 = "lon"
CLASS&lat = lat
CLASS&lon = lon
CLASS@_FillValue       = -1.e+33
replace_ieeenan (CLASS, CLASS@_FillValue, 0)

DOLCE!0 = "lat"
DOLCE!1 = "lon"
DOLCE&lat = lat
DOLCE&lon = lon
DOLCE@_FillValue       = -1.e+33
replace_ieeenan (DOLCE, DOLCE@_FillValue, 0)

CABLE   = CABLE/year_sum
;CABLE(9:10,:,:) = CABLE(9:10,:,:) -273.16
CLASS   = CLASS/year_sum
DOLCE   = DOLCE/year_sum

; ___________________ CERES Net Radiation _________________________
filename3 = case_name(3)+"/CERES_EBAF-Surface_Ed2.8_Subset_200003-201502.nc"
f3 = addfile (filename3,"r")
CERES = f3->sfc_net_tot_all_mon(0,:,:)
;printVarSummary(CERES)
CERES = dim_avg_n_Wrap(lonFlip(f3->sfc_net_tot_all_mon(34:117,:,:)),0)*1.0
; sfc_net_tot_all_mon "Surface Net Total Flux, Monthly Means, All-Sky conditions"
; sfc_net_tot_clr_mon "Surface Net Total Flux, Monthly Means, Clear-Sky conditions"
;printVarSummary(CERES)
;CERES!0 = "lat"
;CERES!1 = "lon"
;CERES&lat = f3->lat
CERES&lon = fspan(0.5,359.5,360)
;printVarSummary(CERES)
; replace_ieeenan (DOLCE, DOLCE@_FillValue, 0)

;_____________________________ PLOT ___________________________________
pic = "spitial_map_vs_obs_off_"+scale
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

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

if scale .eq. "Global" then
  A = -60.
  B = 90.
  C = 0.
  D = 360.
else if scale .eq. "50S-N" then
  A = -50.
  B = 50.
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
end if
res@trYMinF   = A
res@trYMaxF   = B
res@mpMinLatF = A
res@mpMaxLatF = B
res@trXMinF   = C
res@trXMaxF   = D
res@mpMinLonF = C
res@mpMaxLonF = D

if scale .eq. "Global" then
   res@mpCenterLonF       = 160
   res@tmYLValues=(/60,30,0,-30,-60/)
   res@tmYLLabels=(/"60~S~o~N~N","30~S~o~N~N","0~S~o~N~","30~S~o~N~S","60~S~o~N~S"/)
   res@tmXBValues=(/240,300,0,60,120,180/)
   res@tmXBLabels=(/"120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/)
else if scale .eq. "50S-N" then
  res@mpCenterLonF       = 160
  res@tmYLValues=(/50,25,0,-25,-50/)
  res@tmYLLabels=(/"50~S~o~N~N","25~S~o~N~N","0~S~o~N~","25~S~o~N~S","50~S~o~N~S"/)
  res@tmXBValues=(/240,300,0,60,120,180/)
  res@tmXBLabels=(/"120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/)
else if scale .eq. "SE-AU" then
   res@gsnAddCyclic = False
   res@tmYLValues=(/-30,-40/)
   res@tmYLLabels=(/"30~S~o~N~S","40~S~o~N~S"/)
   res@tmXBValues=(/140,150/)
   res@tmXBLabels=(/"140~S~o~N~E","150~S~o~N~E"/)
else if scale .eq. "AU" then
  res@gsnAddCyclic       = False
  res@tmYLValues=(/-10,-20,-30,-40/)
  res@tmYLLabels=(/"10~S~o~N~S","20~S~o~N~S","30~S~o~N~S","40~S~o~N~S"/)
  res@tmXBValues=(/120,130,140,150/)
  res@tmXBLabels=(/"120~S~o~N~E","130~S~o~N~E","140~S~o~N~E","150~S~o~N~E"/)
end if
end if
end if
end if

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
;pres@lbAutoManage       = False
;pres@lbTitleDirection   = "Across"             ; title direction
;pres@lbTitlePosition    = "Bottom";Right"
;pres@lbLabelBarOn       = True


; _________________________ CHECK more _____________________________
filename2 = case_name(1)+"/class_sh_2003-2009.nc"
filename3 = case_name(1)+"/H_estimates_2003-2009.nc"
f1 = addfile (filename2,"r")
f2 = addfile (filename3,"r")

CLASS_sh =  lonFlip(f1->hfss(0,:,:))*1.
CLASS_sh!0 = "lat"
CLASS_sh!1 = "lon"
CLASS_sh&lat = lat
CLASS_sh&lon = lon
CLASS_sh@_FillValue  = -1.e+33
replace_ieeenan (CLASS_sh, CLASS_sh@_FillValue, 0)

CLASS_Pre =  lonFlip(f2->hfss_CLASS_Pre(0,:,:))*1.
CLASS_Pre!0 = "lat"
CLASS_Pre!1 = "lon"
CLASS_Pre&lat = lat
CLASS_Pre&lon = lon
CLASS_Pre@_FillValue  = -1.e+33
replace_ieeenan (CLASS_Pre, CLASS_Pre@_FillValue, 0)

MPI =  lonFlip(f2->hfss_MPI(0,:,:))*1.
MPI!0 = "lat"
MPI!1 = "lon"
MPI&lat = lat
MPI&lon = lon
MPI@_FillValue  = -1.e+33
replace_ieeenan (MPI, MPI@_FillValue, 0)

MERRAFLX =  lonFlip(f2->hfss_MERRAFLX(0,:,:))*1.
MERRAFLX!0 = "lat"
MERRAFLX!1 = "lon"
MERRAFLX&lat = lat
MERRAFLX&lon = lon
MERRAFLX@_FillValue  = -1.e+33
replace_ieeenan (MERRAFLX, MERRAFLX@_FillValue, 0)

MERRALND =  lonFlip(f2->hfss_MERRALND(0,:,:))*1.
MERRALND!0 = "lat"
MERRALND!1 = "lon"
MERRALND&lat = lat
MERRALND&lon = lon
MERRALND@_FillValue  = -1.e+33
replace_ieeenan (MERRALND, MERRALND@_FillValue, 0)

GLDAS =  lonFlip(f2->hfss_GLDAS(0,:,:))*1.
GLDAS!0 = "lat"
GLDAS!1 = "lon"
GLDAS&lat = lat
GLDAS&lon = lon
GLDAS@_FillValue  = -1.e+33
replace_ieeenan (GLDAS, GLDAS@_FillValue, 0)

Princeton =  lonFlip(f2->hfss_Princeton(0,:,:))*1.
Princeton!0 = "lat"
Princeton!1 = "lon"
Princeton&lat = lat
Princeton&lon = lon
Princeton@_FillValue  = -1.e+33
replace_ieeenan (Princeton, Princeton@_FillValue, 0)

a_cor = pattern_cor(CABLE(5,:,:), DOLCE, 1.0, 0) ;unweighted, centered
b_cor = pattern_cor(CABLE(6,:,:), CLASS_Pre, 1.0, 0) ;unweighted, centered
;c_cor = pattern_cor(CABLE(8,:,:), CERES, 1.0, 0) ;unweighted, centered
print("latent heat cor: "+a_cor)
print("sensible heat cor: "+b_cor)
;print("net radiation cor: "+c_cor)


plots = new(8,graphic)

; "Qle" "W/m^2"
  res@cnFillPalette = "WhiteBlueGreenYellowRed";"ViBlGrWhYeOrRe"
  levels = fspan(0,200,21)
  res@cnLevels = levels
  ;pres@cnLevels = levels
  res@tiMainString = "(a) CABLE Latent Heat (W/m~S~2~N~)"
  plots(0) = gsn_csm_contour_map(wks,CABLE(5,:,:),res)
  res@tiMainString = "(b) CLASS Latent Heat (W/m~S~2~N~)"
  plots(1) = gsn_csm_contour_map(wks,DOLCE,res) ;"(f) DOLCE Transpiration (mm/year)"
  ;plots(1) = gsn_csm_contour_map(wks,CLASS(0,:,:),res) ;"(b) CLASS Latent Heat (W/m~S~2~N~)"
  delete(levels)
  delete(res@cnLevels)
  ;delete(pres@cnLevels)

; "Qh" ;"W/m^2"
  levels = fspan(0,200,21)
  res@cnLevels = levels
  ;pres@cnLevels = levels
  res@tiMainString = "(c) CABLE Sensible Heat (W/m~S~2~N~)"
  plots(2) = gsn_csm_contour_map(wks,CABLE(6,:,:),res)
  res@tiMainString = "(d) CLASS Sensible Heat (W/m~S~2~N~)"
  plots(3) = gsn_csm_contour_map(wks,CLASS_Pre,res)
  delete(levels)
  delete(res@cnLevels)
  ;delete(pres@cnLevels)

; "Net Radiation"
  levels = fspan(0,300,31)
  res@cnLevels = levels
  ;pres@cnLevels = levels
;  tmp = CABLE(5,:,:) + CABLE(6,:,:) + CABLE(7,:,:)
;  tmp!0 = "lat"
;  tmp!1 = "lon"
;  tmp&lat = lat
;  tmp&lon = lon
  res@tiMainString = "(e) CABLE Net Radiation (W/m~S~2~N~)"
  plots(4) = gsn_csm_contour_map(wks,CABLE(8,:,:),res)
;  plots(4) = gsn_csm_contour_map(wks,tmp,res)
  ;delete(tmp)
  res@gsnLeftString      = ""
  res@gsnRightString     = ""
  res@tiMainString = "(f) CERES Net Radiation (W/m~S~2~N~)"
  plots(5) = gsn_csm_contour_map(wks,CERES,res)
  delete(levels)
  delete(res@cnLevels)
  ;delete(pres@cnLevels)

;  res@tiMainString = "(f) CABLE  (mm/year)"
;  plots(5) = gsn_csm_contour_map(wks,CABLE(1,:,:),res)

  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = fspan(-100,100,21)
  res@cnLevels = levels
  res@tiMainString = "(g) CABLE Ground Heat (W/m~S~2~N~)"
  plots(6) = gsn_csm_contour_map(wks,CABLE(7,:,:),res)
  delete(levels)
  delete(res@cnLevels)

  res@cnFillPalette = "WhiteBlueGreenYellowRed";
  levels = fspan(0,1.,51)
  res@cnLevels = levels
  res@tiMainString = "(f) Albedo"
  plots(7) = gsn_csm_contour_map(wks,CABLE(9,:,:),res)
  delete(levels)
  delete(res@cnLevels)
  gsn_panel(wks,(/plots/),(/4,2/),pres)

  delete(plots)

  plots = new(6,graphic)

; "Qh" ;"W/m^2"
  res@cnFillPalette = "WhiteBlueGreenYellowRed";
  levels = fspan(0,200,21)
  res@cnLevels = levels
  res@tiMainString = "CLASS_Pre Sensible Heat (W/m~S~2~N~)"
  plots(0) = gsn_csm_contour_map(wks,CLASS_Pre,res)
  res@tiMainString = "MPI Sensible Heat (W/m~S~2~N~)"
  plots(1) = gsn_csm_contour_map(wks,MPI,res)
  res@tiMainString = "MERRAFLX Sensible Heat (W/m~S~2~N~)"
  plots(2) = gsn_csm_contour_map(wks,MERRAFLX,res)
  res@tiMainString = "MERRALND Sensible Heat (W/m~S~2~N~)"
  plots(3) = gsn_csm_contour_map(wks,MERRALND,res)
  res@tiMainString = "GLDAS Sensible Heat (W/m~S~2~N~)"
  plots(4) = gsn_csm_contour_map(wks,GLDAS,res)
  res@tiMainString = "Princeton Sensible Heat (W/m~S~2~N~)"
  plots(5) = gsn_csm_contour_map(wks,Princeton,res)

  gsn_panel(wks,(/plots/),(/3,2/),pres)


; "Qh" ;"W/m^2"
  res@cnFillPalette = "WhiteBlueGreenYellowRed";
  levels = fspan(0,200,21)
  res@cnLevels = levels
  ;pres@cnLevels = levels
  res@tiMainString = "CABLE_off Sensible Heat (W/m~S~2~N~)"
  plots(0) = gsn_csm_contour_map(wks,CABLE(6,:,:),res)
  res@tiMainString = "CLASS_old Sensible Heat (W/m~S~2~N~)"
  plots(1) = gsn_csm_contour_map(wks,CLASS(1,:,:),res)
  res@tiMainString = "CLASS_new Sensible Heat (W/m~S~2~N~)"
  ;tmp = where(CLASS_sh .gt. 0,90,CLASS_sh@_FillValue)
  ;tmp!0 = "lat"
  ;tmp!1 = "lon"
  ;tmp&lat = lat
  ;tmp&lon = lon
  plots(2) = gsn_csm_contour_map(wks,CLASS_sh,res)

  delete(levels)
  delete(res@cnLevels)
  ;delete(pres@cnLevels)

  res@cnFillPalette = "ViBlGrWhYeOrRe";"ViBlGrWhYeOrRe"
  levels = fspan(-60,60,31)
  res@cnLevels = levels
  ;pres@cnLevels = levels

  res@tiMainString = "CLASS_old - CABLE_off"
  tmp = CLASS(1,:,:) - CABLE(6,:,:)
  tmp!0 = "lat"
  tmp!1 = "lon"
  tmp&lat = lat
  tmp&lon = lon
  printVarSummary(tmp)
  plots(3) = gsn_csm_contour_map(wks,tmp,res)
  delete(tmp)

  res@tiMainString = "CLASS_new - CABLE_off"
  tmp = CLASS_sh - CABLE(6,:,:)
  tmp!0 = "lat"
  tmp!1 = "lon"
  tmp&lat = lat
  tmp&lon = lon
  printVarSummary(tmp)
  plots(4) = gsn_csm_contour_map(wks,tmp,res)
  delete(tmp)

  res@tiMainString = "CLASS_new - CLASS_old"
  tmp = CLASS_sh - CLASS(1,:,:)
  tmp!0 = "lat"
  tmp!1 = "lon"
  tmp&lat = lat
  tmp&lon = lon
  printVarSummary(tmp)
  plots(5) = gsn_csm_contour_map(wks,tmp,res)
  delete(tmp)
  delete(levels)
  delete(res@cnLevels)
  ;delete(pres@cnLevels)

  gsn_panel(wks,(/plots/),(/3,2/),pres)


end
