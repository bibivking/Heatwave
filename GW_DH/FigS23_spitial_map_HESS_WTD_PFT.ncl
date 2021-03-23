
;MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
; Author  : Mengyuan Mu
; Version : 1.0 (22.03.2021)"
; Email   : mu.mengyuan815@gmail.com
;WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

;*******************************************
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
year_s      = 2017
year_e      = 2019
scale       = "SE-AU"

PFT         = new((/241,281/),integer,-9999999)
WTD         = new((/241,281/),float)
PFT         = -9999999
WTD         = 0.

do year  = year_s,year_e
  filename_hw = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/CTL-30x3+30yr/outputs-daily/cable_out_"+year+"_SE_Aus.nc"
  f           = addfile (filename_hw,"r")
  WTD         = WTD + dim_avg_n_Wrap(f->WatTable(:,:,:),0)
  ; delete(f)
end do
PFT = f->iveg
PFT = where(PFT .eq. 2, 1, PFT)
PFT = where(PFT .eq. 5, 2, PFT)
PFT = where(PFT .eq. 6, 3, PFT)
PFT = where(PFT .eq. 9, 4, PFT)
PFT = where(PFT .eq. 14,5, PFT)
print(get_unique_values(PFT))

WTD = WTD/3.
WTD@_FillValue = f->WatTable@_FillValue
WTD!0 = "lat"
WTD!1 = "lon"
WTD&lat = f->y ;fspan(-40,-28,241)
WTD&lon = f->x ;fspan(140,154,281)
printVarSummary(WTD)

; ================================= Plotting ===================================
res                    = True
res@gsnDraw            = True  ; Don't draw plots
res@gsnFrame           = True  ; ����ҳ
res@gsnMaximize        = False

res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                    ;��ɫ
res@lbLabelBarOn           = True
res@vpHeightF 	           = 0.38               ; change aspect ratio of plot
res@vpWidthF 	             = 0.38
res@txFontHeightF          = 0.015            ; change title font heights
res@tiYAxisFontThicknessF  = 3.
res@tmBorderThicknessF     = 3.
res@tmXBMajorThicknessF    = 3.
res@tmXBMinorThicknessF    = 3.
res@tmYLMajorThicknessF    = 3.
res@tmYLMinorThicknessF    = 3.
res@tmYLLabelFontThicknessF= 3.
res@tmBorderThicknessF     = 3.


;************** ����labelbar ***************
res@lbBoxLinesOn       = True                       ;�ر�lbar box ����
; res@lbTitleFont        = 25
; res@lbLabelFont        = 25
res@lbTitleFontHeightF = 0.015
res@lbLabelFontHeightF = 0.015
res@txString           = ""
; res@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
; res@tmYLLabelFont      = 25
res@pmLabelBarHeightF       = 0.08
res@pmLabelBarParallelPosF  = 0.5
res@pmLabelBarOrthogonalPosF= 0.09

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
res@tmXBLabelFontHeightF    = 0.015
res@tmYLLabelFontHeightF    = 0.015
res@tmYLMode                = "Explicit"
res@tmXBMode                = "Explicit"

res@gsnMaximize         = True
res@gsnLeftString       = ""
res@gsnRightString      = ""


res@cnFillMode             = "RasterFill"            ; Raster Mode
res@cnFillOn               = True                            ;��ɫ
res@cnLinesOn              = False
res@cnLineLabelsOn         = False
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode   = "ExplicitLevels"

res@lbLabelBarOn           = True
res@lbBoxLinesOn           = True                       ;�ر�lbar box ����
res@lbOrientation          = "horizontal";"vertical" ;"horizontal"
; res@lbTitleFont            = 25
; res@lbLabelFont            = 25
res@lbTitleFontHeightF     = 0.013
res@lbLabelFontHeightF     = 0.013

; res@tmXBLabelFont           = 25 ;Sets the font index for the bottom X-Axis labels.
; res@tmYLLabelFont           = 25
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF    = 0.015
res@tmYLLabelFontHeightF    = 0.015
; res@mpDataBaseVersion       = "HighRes"

res@lbLabelBarOn            = True
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
  res@tmYLLabels      = (/"40~S~o~N~S","35~S~o~N~S","30~S~o~N~S","25~S~o~N~S","20~S~o~N~S","15~S~o~N~S","10~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/115,120,125,130,135,140,145,150/)
  res@tmXBLabels      = (/"115~S~o~N~E","120~S~o~N~E","125~S~o~N~E","130~S~o~N~E","135~S~o~N~E","140~S~o~N~E","145~S~o~N~E","150~S~o~N~E"/)
end if
end if
end if

pic1 = "./plots/WTD_2017-2019"
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"GMT_haxby");"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

levels = fspan(1,20,20)
res@cnLevels = levels
; res@txString = "Water Table Depth (CTL-30+30yr)" ;var@long_name
plots = gsn_csm_contour_map(wks1,WTD,res)
delete(levels)
delete(res@cnLevels)

pic2 = "./plots/PFT"
wks2 = gsn_open_wks("pdf",pic2)
gsn_define_colormap(wks2,"drought_severity");"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")
levels = (/1.1,2.1,3.1,4.1/)
res@lbLabelsOn                 = True
res@cnExplicitLabelBarLabelsOn = True
res@lbLabelStrings             = (/"evergreen~C~broadleaf~C~forest",\
                                   "shrub~C~ ~C~ ", "grassland~C~ ~C~ ", "cropland~C~ ~C~ ",\
                                   "barren~C~land~C~ "/)
res@lbTitleFontHeightF         = 0.01
res@lbLabelFontHeightF         = 0.01
res@pmLabelBarParallelPosF     = 0.5
res@pmLabelBarOrthogonalPosF   = 0.13
; res@lbLabelPosition            = "right"
res@lbLabelJust                = "bottomcenter"
res@lbLabelAlignment           = "BoxCenters"
; res@lbLabelOffsetF             = -0.1
res@cnLevels = levels
plots = gsn_csm_contour_map(wks2,PFT,res)

end
