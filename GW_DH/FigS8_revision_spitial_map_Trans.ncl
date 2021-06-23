
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

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 500000000000 ;
end setvalues
;_____________________________________

; =========================== Notes for time step data =========================
; model time step : 0, 1, 2, 3 , 4 , 5,  6,  7
; real world time : 0, 3, 6, 9, 12, 15, 18, 21
; ==============================================================================

begin

; =========================== OPTIONS ===========================
PFT        = "tree" ; "tree"
scale      = "SE-AU"
year_s     = 2019
year_e     = 2019

HW_day     = (/ 15, 25 /)
event_sum  = dimsizes(HW_day)

message    = "Trans"

pic_name  = message+"_"+scale+"_"+year_s+"-"+year_e

path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"

case_name = (/"CTL-30x3+30yr","FREE_DRAIN","CTL-30x3+30yr_2019-DeepRoot"/)
case_sum  = dimsizes(case_name)
case_name_brf = (/"GW","FD","DR"/)

; =========================== Pre-define ===========================
if PFT .eq. "all" then
  iveg_s = 0
  iveg_e = 20
else if PFT .eq. "tree" then
  iveg_s = 1
  iveg_e = 4
end if
end if

; region: AU
if scale .eq. "AU" then
  A = -44.
  B = -10.
  C = 112.
  D = 154.
else if scale .eq. "SE-AU" then
  A = -40.
  B = -28.
  C = 140.
  D = 154.
end if
end if

; sub region 1
A_rgn_1 = -36
B_rgn_1 = -34
C_rgn_1 = 149
D_rgn_1 = 151

; sub region 2
A_rgn_2 = -38
B_rgn_2 = -36
C_rgn_2 = 146.5
D_rgn_2 = 148.5

; ______ load path ______
path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = path_mum+case_name(case_num)+"/outputs-subday"
end do

filename1 = path(0)+"/cable_out_2019_SE_Aus.nc"
f1        = addfile (filename1,"r")
lat       = f1->y({A:B})
lon       = f1->x({C:D})
latitude  = f1->latitude({A:B},{C:D})
longitude = f1->longitude({A:B},{C:D})
lat_sum   = dimsizes(lat)
lon_sum   = dimsizes(lon)
default_val = f1->Qle@_FillValue
delete(f1)
delete(filename1)

order_name  = (/"(a)","(b)","(c)","(d)",  "(e)","(f)","(g)","(h)",\
                "(i)","(j)","(k)","(l)",  "(m)","(n)","(o)","(p)",\
                "(q)","(r)","(s)","(t)"/)

; ================================== TVeg ======================================
TVeg     = new((/case_sum,2,lat_sum,lon_sum/),float,default_val) ; case x days x lat x lat
TVeg     = 0.

do case_num = 0, case_sum -1
  print("case_num = " +case_num)
  do year = year_s,year_e
    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
    f1 = addfile (filename1,"r")

    TVeg_tmp   = f1->TVeg(0:60*8,{A:B},{C:D})
    do event_num = 0,event_sum-1
      ; 0 am
      ts = (HW_day(event_num)-1)*8
      ; to compare with Aqua  - 2pm
      TVeg(case_num, event_num, :, :) = TVeg_tmp(ts+4,:,:) + (TVeg_tmp(ts+5,:,:) - TVeg_tmp(ts+4,:,:))*(2/3)
    end do
  end do ; year
end do ; case_num

TVeg        = TVeg*3600.*24.*28.94 ; mm/s -> mm/d -> W/m2, 1 mm/day = 28.94 W/m2
TVeg!0      = "case"
TVeg!1      = "time"
TVeg!2      = "lat"
TVeg!3      = "lon"
TVeg&lat    = lat
TVeg&lon    = lon

; ================================ Calc Diff ===================================
TVeg_GW_FD = TVeg(0,:,:,:)
TVeg_GW_FD = TVeg_GW_FD - TVeg(1,:,:,:)
TVeg_DR_GW = TVeg(2,:,:,:)
TVeg_DR_GW = TVeg_DR_GW - TVeg(0,:,:,:)


; ================================= Plotting ===================================
pic = "./plots/FigS8_revision_"+pic_name
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")

res                         = True
res@gsnDraw                 = False  ; Don't draw plots
res@gsnFrame                = False  ; ����ҳ
res@gsnMaximize             = False

res@cnFillOn                = True
res@cnMissingValFillColor   = "gray80"
res@cnFillBackgroundColor   = "gray80"

res@vpHeightF 	            = 0.38               ; change aspect ratio of plot
res@vpWidthF 	              = 0.38
res@tiYAxisFontThicknessF   = 3.
res@tmXBMajorThicknessF     = 3.
res@tmXBMinorThicknessF     = 3.
res@tmYLMajorThicknessF     = 3.
res@tmYLMinorThicknessF     = 3.
res@tmYLLabelFontThicknessF = 3.
res@tmBorderThicknessF      = 3.
; res@txFont                  = 25
res@txFontHeightF           = 0.015            ; change title font heights
res@gsnPaperOrientation     = "landscape";"portrait";"landscape" ;"auto"

res@cnFillMode             = "RasterFill"            ; Raster Mode
res@cnFillOn               = True                            ;��ɫ
res@cnLinesOn              = False
res@cnLineLabelsOn         = False
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode   = "ExplicitLevels"

res@lbLabelBarOn           = True
res@lbBoxLinesOn           = True                       ;�ر�lbar box ����
res@lbOrientation          = "vertical";"vertical" ;"horizontal"
res@lbTitleFont            = 25
res@lbLabelFont            = 25
res@lbTitleFontHeightF     = 0.02
res@lbLabelFontHeightF     = 0.015

res@tmXBLabelFont           = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont           = 25
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF    = 0.015
res@tmYLLabelFontHeightF    = 0.015
; res@mpDataBaseVersion       = "HighRes"

res@lbLabelBarOn            = True
res@lbTitlePosition         = "Right"
res@lbTitleAngleF           = 270
res@lbTitleDirection        = "Across"                          ; letter angle
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

; _____________ Colorbar ____________
; Red
cmap1  = (/ (/255.,255,255,1/), (/255,255,229,1/), (/255,247,188,1/), (/254,227,145,1/), \
            (/254,196,79,1/),  (/254,153,41,1/),  (/236,112,20,1/),  (/204,76,2,1/),    \
            (/153,52,4,1/),    (/139,35,35,1/), (/102,37,6,1/)    /)
; Blue
cmap2  = (/ (/255.,255,255,1/), (/247,252,240,1/), (/224,243,219,1/), (/204,235,197,1/),  \
            (/168,221,181,1/), (/123,204,196,1/),  (/78,179,211,1/),  (/43,140,190,1/),  \
            (/8,104,172,1/),    (/  8,64,129,1/),    (/2,56,88,1/) /)

; Blue to Red
cmap3  = (/ (/8.,29,88,1/),    (/  8,64,129,1/),  (/8,104,172,1/),  \                   \
            (/43,140,190,1/),  (/78,179,211,1/),  (/123,204,196,1/),(/168,221,181,1/), \
            (/204,235,197,1/), (/224,243,219,1/), (/255,255,255,1/),\
            (/255,255,255,1/), (/255,247,188,1/), (/254,227,145,1/),\
            (/254,196,79,1/),  (/254,153,41,1/),  (/236,112,20,1/),  (/204,76,2,1/),   \
            (/153,52,4,1/),    (/139,35,35,1/),    (/102,37,6,1/)   /) ;(/128,0,38,1/)   /) ;  (/102,37,6,1/) (/103,0,13,1/) (/2.,56,88,1/),


cmap1(:,0:2) = cmap1(:,0:2)/255.
cmap2(:,0:2) = cmap2(:,0:2)/255.
cmap3(:,0:2) = cmap3(:,0:2)/255.

;**************** ͼ�������� *****************
pres                    = True
pres@wkWidth            = 2500
pres@wkHeight           = 2500
pres@vpHeightF 	        = 0.2               ; change aspect ratio of plot
pres@vpWidthF 	        = 0.95
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False
; pres@gsnPanelMainPosXF  = 0.1
; pres@gsnPanelMainPosYF  = 0.1
pres@cnLevelSelectionMode = "ExplicitLevels"
pres@lbBoxLinesOn       = True
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbTitleFontHeightF = 0.02
pres@lbLabelFontHeightF = 0.02
pres@gsnPanelYWhiteSpacePercent = 0.5
pres@gsnPanelXWhiteSpacePercent = 0.5


; -------------- plot order ----------------
txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.03

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

; -------------- start to plot ----------------
plots  = new(10,graphic)

; -------------- FigS8 a-b GW-FD --------------
res@lbTitleString = "~F25~Et~B~GW_2pm~N~ - ~F25~Et~B~FD_2pm~N~ (W m~S~-2~N~)"
res@cnFillPalette = cmap3(::-1,:) ;"precip3_16lev";"testcmap";"ViBlGrWhYeOrRe";"WhiteBlueGreenYellowRed"
res@cnLevels      = (/-250,-200,-150,-100,-50,-25,25,50,100,150,200,250/)
res@lbLabelBarOn = False ; "Always"
plots(0) = gsn_csm_contour_map(wks,TVeg_GW_FD(0,:,:),res) ;LST_Aqua_regrid(14,:,:),res)

res@lbLabelBarOn = True
plots(1) = gsn_csm_contour_map(wks,TVeg_GW_FD(1,:,:),res)


; -------------- FigS6 c-d DR-GW --------------
; notes :
; VegT: case x 2 time steps x days x lat x lon
res@lbTitleString =  "~F25~Et~B~DR_2pm~N~ - ~F25~Et~B~GW_2pm~N~ (W m~S~-2~N~)"
res@lbLabelBarOn = False
plots(2) = gsn_csm_contour_map(wks,TVeg_DR_GW(0,:,:),res)
res@lbLabelBarOn = True
plots(3) = gsn_csm_contour_map(wks,TVeg_DR_GW(1,:,:),res)
delete(res@cnLevels)
delete(res@cnFillPalette)

; ________________ add shading ________________
poly             = new(4,graphic)
gres             = True
gres@gsEdgesOn   = True
gres@gsEdgeColor = "red"
gres@gsLineColor = "red"
gres@gsLineThicknessF = 3.0

xmin = (/C_rgn_1,C_rgn_2/)
xmax = (/D_rgn_1,D_rgn_2/)
ymin = (/A_rgn_1,A_rgn_2/)
ymax = (/B_rgn_1,B_rgn_2/)

xbars = (/ (/xmin(0), xmax(0), xmax(0), xmin(0), xmin(0)/), \
            (/xmin(1), xmax(1), xmax(1), xmin(1), xmin(1)/) /)

ybars = (/ (/ymin(0),ymin(0),ymax(0),ymax(0),ymin(0)/), \
           (/ymin(1),ymin(1),ymax(1),ymax(1),ymin(1)/)  /)

poly(0) = gsn_add_polyline(wks,plots(0),xbars(0,:),ybars(0,:),gres)
poly(1) = gsn_add_polyline(wks,plots(1),xbars(1,:),ybars(1,:),gres)

poly(2) = gsn_add_polyline(wks,plots(2),xbars(0,:),ybars(0,:),gres)
poly(3) = gsn_add_polyline(wks,plots(3),xbars(1,:),ybars(1,:),gres)

do cnt = 0, 3
  txid1   = gsn_create_text(wks, order_name(cnt), txres)
  annoid1 = gsn_add_annotation(plots(cnt), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/5,2/),pres)

end
