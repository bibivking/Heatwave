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

; ================================== Options ===================================
year_s        = 2017
year_e        = 2019
year_sum      = year_e - year_s+1
time_tot      = 365 ; *8
scale         = "SE-AU"

; filename_hw = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/CTL-30x3+30yr/outputs-daily/cable_out_"+year_s+"_SE_Aus.nc"
; f           = addfile (filename_hw,"r")
; ; var1  = dim_sum_n_Wrap(f->Qrecharge(:,:,:)*24*3600,0)/year_sum
; iveg   = f->iveg(:,:)
; var1   = where( iveg .ge. 1 .and. iveg .le. 4, dim_avg_n_Wrap(f->LAI(:,:,:),0), f->LAI@_FillValue)

WTD       = new((/241,281/),float)
var1      = new((/241,281/),float)
var2      = new((/241,281/),float)
WTD       = 0.
var1      = 0.
var2      = 0.

do year = year_s, year_e
  print("year="+year)
  filename_gw = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/CTL-30x3+30yr/outputs-daily/cable_out_"+year+"_SE_Aus.nc"
  f           = addfile (filename_gw,"r")
  default_val = f->SoilMoist@_FillValue
  var1        = var1 + dim_avg_n_Wrap(f->SoilMoist(:,5,:,:),0)/year_sum ;*24*2600
  WTD         = WTD + dim_avg_n_Wrap(f->WatTable(:,:,:),0)/year_sum
  ; delete(f)
  filename_fd = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/FREE_DRAIN/outputs-spinup30x3yr/cable_out_"+year+"_SE_Aus.nc"
  f           = addfile (filename_fd,"r")
  var2        = var2 + dim_avg_n_Wrap(f->SoilMoist(:,5,:,:),0)/year_sum ;*24*2600
  ; delete(f)
end do

var1@_FillValue = default_val
var2@_FillValue = default_val

deltaSM = f->SoilMoist(0,5,:,:)
deltaSM = var2*1. ;var1 - var2

WTD@_FillValue = default_val
WTD!0 = "lat"
WTD!1 = "lon"
WTD&lat = f->y ;fspan(-40,-28,241)
WTD&lon = f->x ;fspan(140,154,281)
printVarSummary(WTD)


; -------------------------------------------
WTD_class     = new((/5, 241,281/),float, default_val)
WTD_class@_FillValue = -9999.
WTD_class!0 = "WTD_lvl"
WTD_class!1 = "lat"
WTD_class!2 = "lon"
WTD_class&lat = f->y ;fspan(-40,-28,241)
WTD_class&lon = f->x ;fspan(140,154,281)

SM6_WTD_class = new((/5, 241,281/),float, default_val)
SM6_WTD_class@_FillValue = -9999.
SM6_WTD_class!0 = "WTD_lvl"
SM6_WTD_class!1 = "lat"
SM6_WTD_class!2 = "lon"
SM6_WTD_class&lat = f->y ;fspan(-40,-28,241)
SM6_WTD_class&lon = f->x ;fspan(140,154,281)
; -------------------------------------------
WTD_class(0,:,:) = where(WTD .le. 4,                   WTD,WTD@_FillValue)
WTD_class(1,:,:) = where(WTD .le. 4.6,                 WTD,WTD@_FillValue)
WTD_class(2,:,:) = where(WTD .gt. 4.6 .and. WTD .lt. 5,WTD,WTD@_FillValue)
WTD_class(3,:,:) = where(WTD .ge. 4.6 .and. WTD .lt. 7, WTD,WTD@_FillValue)
WTD_class(4,:,:) = where(WTD .ge. 10,                  WTD,WTD@_FillValue)


SM6_WTD_class(0,:,:) = where(WTD .le. 4,                   deltaSM,deltaSM@_FillValue)
SM6_WTD_class(1,:,:) = where(WTD .le. 4.6,                 deltaSM,deltaSM@_FillValue)
SM6_WTD_class(2,:,:) = where(WTD .gt. 4.6 .and. WTD .lt. 5,deltaSM,deltaSM@_FillValue)
SM6_WTD_class(3,:,:) = where(WTD .ge. 4.6 .and. WTD .lt. 7,deltaSM,deltaSM@_FillValue)
SM6_WTD_class(4,:,:) = where(WTD .ge. 10,                  deltaSM,deltaSM@_FillValue)



; filename_hw = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/FREE_DRAIN/outputs-spinup30x3yr/cable_out_"+year_s+"_SE_Aus.nc"
; f           = addfile (filename_hw,"r")
; var  = dim_sum_n_Wrap(f->Qrecharge(:,:,:)*24*3600,0)/year_sum
;
; do year = year_s+1, year_e
; filename_hw = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/FREE_DRAIN/outputs-spinup30x3yr/cable_out_"+year+"_SE_Aus.nc"
; f           = addfile (filename_hw,"r")
; var         = var + dim_sum_n_Wrap(f->Qrecharge(:,:,:)*24*3600,0) ;*24*2600
; delete(f)
; end do

;
; filename_hw = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/CTL-30x3+30yr-DeepRoot/outputs-daily/cable_out_"+year+"_SE_Aus.nc"
; f           = addfile (filename_hw,"r")
; var2 = f->froot(:,:,:)
; delete(f)

; var2 = new((/12,12,14/),float)
; do month = 1,12
;   if month .lt. 10 then
;     filename_hw = "/g/data/w35/mm3972/data/CERES/CERES_EBAF_Ed4.1_Subset_20130"+month+"-"+"20130"+month+".nc"
;   else
;     filename_hw = "/g/data/w35/mm3972/data/CERES/CERES_EBAF_Ed4.1_Subset_2013"+month+"-"+"2013"+month+".nc"
;   end if
;
;   f           = addfile (filename_hw,"r")
;   print(f->sfc_net_tot_clr_t_mon(0,{-40:-28},{140:154}))
;   var2(month-1,:,:) = f->sfc_net_tot_clr_t_mon(0,{-40:-28},{140:154})
;   delete(f)
; end do

; filename_hw = "/g/data/w35/Shared_data/AWAP_3h_v1/SWdown/AWAP.SWdown.3hr.2013.nc"
; f           = addfile (filename_hw,"r")
; var2 = f->SWdown(:,{-40:-28},{140:154})
; delete(f)


; filename_hw = "/g/data/w35/Shared_data/Observations/AWAP_all_variables/daily/rad/AWAP_daily_rad_1990_2019.nc"
; f           = addfile (filename_hw,"r")
; var2        = f->rad(:,{-40:-28},{140:154})
; delete(f)
; ================================= Plotting ===================================
res                         = True
res@gsnDraw                 = True  ; Don't draw plots
res@gsnFrame                = True  ; ����ҳ
res@gsnMaximize             = True

res@vpHeightF 	            = 0.4               ; change aspect ratio of plot
res@vpWidthF 	              = 0.8
res@tiYAxisFontThicknessF   = 3.
res@tmXBMajorThicknessF     = 3.
res@tmXBMinorThicknessF     = 3.
res@tmYLMajorThicknessF     = 3.
res@tmYLMinorThicknessF     = 3.
res@tmYLLabelFontThicknessF = 3.
res@tmBorderThicknessF      = 3.
res@txFont                  = 25
res@txFontHeightF           = 0.015            ; change title font heights
res@gsnPaperOrientation     = "landscape";"portrait";"landscape" ;"auto"

res@cnFillMode             = "RasterFill"            ; Raster Mode
res@cnFillOn               = True                            ;��ɫ
res@cnLinesOn              = False
res@cnLineLabelsOn         = False
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode   = "ExplicitLevels"

res@lbLabelBarOn            = True
res@lbBoxLinesOn            = True                       ;�ر�lbar box ����
res@lbTitleFont             = 25
res@lbLabelFont             = 25
res@lbTitleFontHeightF      = 0.013
res@lbLabelFontHeightF      = 0.013
; res@lbOrientation           = "Vertical"
res@lbTitlePosition         = "Bottom"
res@pmLabelBarHeightF       = 0.1
res@pmLabelBarParallelPosF  = 0.5
res@pmLabelBarOrthogonalPosF= 0.09

res@tmXBLabelFont           = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont           = 25
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF    = 0.015
res@tmYLLabelFontHeightF    = 0.015
; res@mpDataBaseVersion       = "HighRes"

res@tmYLMode                = "Explicit"
res@tmXBMode                = "Explicit"
res@txString                = ""
res@gsnLeftString           = ""
res@gsnRightString          = ""

if scale .eq. "Global" then
  res@trYMinF         = -60
  res@trYMaxF         = 80
  res@mpMinLatF       = -60
  res@mpMaxLatF       = 80
  res@mpCenterLonF    = 160
  res@tmYLMode        = "Explicit"
  res@tmYLValues      =(/60,30,0,-30,-60/)
  res@tmYLLabels      =(/"60~S~o~N~N","30~S~o~N~N","0~S~o~N~","30~S~o~N~S","60~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      =(/240,300,0,60,120,180/)
  res@tmXBLabels      =(/"120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/)
else if scale .eq. "SE-AU" then
  res@gsnAddCyclic    = False
  res@trYMinF         = -40.
  res@trYMaxF         = -28.
  res@mpMinLatF       = -40.
  res@mpMaxLatF       = -28.
  res@trXMinF         = 140.
  res@trXMaxF         = 154.
  res@mpMinLonF       = 140.
  res@mpMaxLonF       = 154.
  res@tmYLMode        = "Explicit"
  res@tmYLValues      = (/-40,-35,-30,-25/)
  res@tmYLLabels      = (/"40~S~o~N~S","35~S~o~N~S","30~S~o~N~S","25~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/140,145,150/)
  res@tmXBLabels      = (/"140~S~o~N~E","145~S~o~N~E","150~S~o~N~E"/)
else if scale .eq. "AU" then
  res@gsnAddCyclic    = False
  res@trYMinF         = -44
  res@trYMaxF         = -10
  res@mpMinLatF       = -44
  res@mpMaxLatF       = -10
  res@trXMinF         = 112
  res@trXMaxF         = 154
  res@mpMinLonF       = 112
  res@mpMaxLonF       = 154
  res@tmYLMode        = "Explicit"
  res@tmYLValues      = (/-40,-35,-30,-25,-20,-15,-10/)
  res@tmYLLabels      = (/"40~S~o~N~N","35~S~o~N~N","30~S~o~N~","25~S~o~N~S","20~S~o~N~N","15~S~o~N~N","10~S~o~N~"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/115,120,125,130,135,140,145,150/)
  res@tmXBLabels      = (/"115~S~o~N~W","120~S~o~N~W","125~S~o~N~","130~S~o~N~W","135~S~o~N~W","140~S~o~N~","145~S~o~N~W","150~S~o~N~"/)
end if
end if
end if

res@gsnLeftString       = ""
res@gsnRightString      = ""

pic1 = "./plots/Check_SM=6_WTD_class_2017-2019"
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"spread_15lev");"spread_15lev");"WhiteBlueGreenYellowRed");"GMT_haxby");"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

res@gsnLeftString = "Soil Moisture difference" ;"Water Table Depth" ;var@long_name
;levels = fspan(-0.2,0.2,11) ;fspan(-500.,500.,51)

levels = fspan(0.,0.5,11) ;fspan(-500.,500.,51)
res@cnLevels = levels
; do soil = 0, 5
plots = gsn_csm_contour_map(wks1,deltaSM,res)

res@gsnLeftString = "WTD < 4"
plot              = gsn_csm_contour_map(wks1,SM6_WTD_class(0,:,:), res)
res@gsnLeftString = "WTD < 4.6"
plot              = gsn_csm_contour_map(wks1,SM6_WTD_class(1,:,:), res)
res@gsnLeftString = "4.6 < WTD < 5"
plot              = gsn_csm_contour_map(wks1,SM6_WTD_class(2,:,:), res)
res@gsnLeftString = "4.6 < WTD < 7"
plot              = gsn_csm_contour_map(wks1,SM6_WTD_class(3,:,:), res)
res@gsnLeftString = "10 < WTD"
plot              = gsn_csm_contour_map(wks1,SM6_WTD_class(4,:,:), res)

delete(levels)
delete(res@cnLevels)
levels = fspan(1,15,15) ;fspan(-500.,500.,51)
res@cnLevels = levels
res@gsnLeftString = "WTD < 4"
plot              = gsn_csm_contour_map(wks1,WTD_class(0,:,:), res)
res@gsnLeftString = "WTD < 4.6"
plot              = gsn_csm_contour_map(wks1,WTD_class(1,:,:), res)
res@gsnLeftString = "4.6 < WTD < 5"
plot              = gsn_csm_contour_map(wks1,WTD_class(2,:,:), res)
res@gsnLeftString = "4.6 < WTD < 7"
plot              = gsn_csm_contour_map(wks1,WTD_class(3,:,:), res)
res@gsnLeftString = "10 < WTD"
plot              = gsn_csm_contour_map(wks1,WTD_class(4,:,:), res)


; end do
  ; plots = gsn_csm_contour(wks1,var2(day_num, :,:),res)
; end do
; counter = 0
; do year = 1990, 2019
;   if mod(year,4) .eq. 0 then
;     day_sum = 366
;   else
;     day_sum = 365
;   end if
;   do day_num = 0, day_sum -1
;     res@gsnLeftString = year+"-"+day_num
;     if day_num .lt. 60 .or. day_num .gt. day_sum - 61 then
;       plots = gsn_csm_contour(wks1,var2(counter,:,:),res)
;       counter = counter + 1
;     else
;       counter = counter + 1
;     end if
;   end do
; end do

; delete(levels)
; delete(res@cnLevels)
; levels = fspan(-10,10,21)
; res@cnLevels = levels
; plots = gsn_csm_contour(wks1,var-var2,res)
end
