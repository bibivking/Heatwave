
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
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

begin

; =============================== Pre-define ===================================
case_name   = (/"CTL","FD"/) ;,"DEEP_ROOT"/)
path        = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
pic_message = "CTL+FD_in_1_panel"

year_s         = 2000
year_e         = 2019

year_sum       = year_e-year_s+1
total          = year_sum*12
case_sum       = dimsizes(case_name)
default_value  = -99999.

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)
order_name     = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)

; ================================ Set Regrid ==================================
PFT           = "all" ; "tree" ;"all"
scale         = "SE-AU"
time_scale    = "monthly" ; "annual"; "monthly"
cable_version = "CABLE-HESS" ;"LIS-CABLE" ; "CABLE-HESS"
grace_version = "Original-fill-nan"
                ;"Original-fill-nan": using Nan fill missing month
                ;"Gap-fill" : Sanaa's gap fill
                ;"Original" : if this option, should alwasy set GRACE_start_year = 2003
                ;             GRACE_end_year <= 2010, since too many gaps after 2010


; ================================= PLOT ======================================
pic = "./plots/Fig1a_time_serial_GRACE_JPL06-"+grace_version+"_vs_"+cable_version+\
      "_PFT-"+PFT+"_"+pic_message+"_"+scale+"_"+time_scale
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")

res                        = True              ; plot mods desired
res@gsnDraw                = False             ; don't draw yet
res@gsnFrame               = False             ; don't advance frame yet
res@gsnMaximize            = False
res@vpHeightF 	           = 0.4               ; change aspect ratio of plot
res@vpWidthF 	             = 0.8
res@txFontHeightF          = 0.015            ; change title font heights
res@tiYAxisFontThicknessF  = 2.
res@tmXBMajorThicknessF    = 2.
res@tmXBMinorThicknessF    = 2.
res@tmYLMajorThicknessF    = 2.
res@tmYLMinorThicknessF    = 2.
res@tmYLLabelFontThicknessF= 2.
res@tmBorderThicknessF     = 2.

; res@txFont                 = "times-bold"
res@xyMonoLineColor        = False             ; want colored lines
res@tmXBMode               = "Explicit"
; res@tiXAxisString          = "Year"
res@tfPolyDrawOrder        = "PreDraw"

res@lgPerimOn              = False      ; No legend perimeter.
; res@lgJustification        = "TopRight"
; res@pmLegendSide           = "Bottom"
res@pmLegendParallelPosF   = 0.7       ; Move to right
res@pmLegendOrthogonalPosF = -0.55     ; Move into plot
res@pmLegendWidthF         = 0.15       ; Decrease width
res@pmLegendHeightF        = 0.11       ; Decrease height
res@lgBoxMinorExtentF      = 0.3       ; Shorten the legend lines

LineColors                 = (/ "black","darkgreen","cornflowerblue","orange"/)
LineThicknesses	           = (/2.5,2.5,1.5,2.5/)      ; line thicknesses
DashPatterns	             = (/0.,0.,0.,0./)      ; make all lines solid
Opacities                  = (/0.9,0.9,0.9,0.9/)
Labels                     = (/"GRACE","GW", "GW Aquifer", "FD"/)
res@lgItemOrder            = (/ 2,3,1,0 /)   ; Reorder the legends

res@xyLineColors           = LineColors
res@xyLineThicknesses      = LineThicknesses
res@xyDashPatterns	       = DashPatterns
res@xyExplicitLegendLabels = Labels
res@xyLineOpacities        = Opacities
x                          = ispan(1,total,1)
res@tmXBValues             = ispan(1,total+1,12*2)
res@tmXBLabels             = tostring(ispan(year_s,year_e+1,2))
res@trXMaxF                = max(x)+1
res@trXMinF                = min(x)
res@trYMaxF                = 200.
res@trYMinF                = -100.


pres                            = True
pres@wkWidth                    = 2500
pres@wkHeight                   = 2500
pres@gsnMaximize                = False
pres@gsnPanelYWhiteSpacePercent = 0.1
pres@gsnPanelXWhiteSpacePercent = 0.1

plots = new(1,graphic)
poly  = new(2,graphic)

res@tiYAxisString          = "TWSA (mm)"
res@pmLegendDisplayMode    = "Always"

; ________________ read ascii ________________
; In the txt file, 0: GRACE; 1: TWS; 2 GWMoist;
;                  3 SoilMoist; 4 CanopInt; 5 SWE
tmp_all           = new((/4,total/),"float") ; 0: GRACE, 1: GW TWS, 2: GW Aquifer 3: FD SM

; GW
txt_name_GW       = "./txt/"+case_name(0)+"_GRACE_PFT-"+PFT+"_"+year_s+"-"+year_e+".txt"
tmp_GW            = asciiread(txt_name_GW,(/6,total/),"float")
tmp_GW@_FillValue = default_value
SD_GRACE          = stddev(tmp_GW(0,:))
r_GW              = escorc(tmp_GW(1,:), tmp_GW(0,:))   ; Pearson’s correlation coefficient (r)

GW_GRACE_with_data= where( .not. ismissing(tmp_GW(0,:)), tmp_GW(1,:), default_value)
GW_GRACE_with_data@_FillValue = default_value
SD_GW             = stddev(GW_GRACE_with_data) ;stddev(tmp_GW(1,:)) ; TWS
MBE_GW            = avg(tmp_GW(1,:) - tmp_GW(0,:))     ; Mean Bias Error (MBE)
tmp_all(0:2,:)    = tmp_GW(0:2,:)
delete(tmp_GW)

; FD, excluding aquifer when calculating TWSA, since aquifer shouldn't have an impact in FD
txt_name_FD       = "./txt/"+case_name(1)+"_GRACE_PFT-"+PFT+"_"+year_s+"-"+year_e+".txt"
tmp_FD            = asciiread(txt_name_FD,(/6,total/),"float")
tmp_FD@_FillValue = default_value
r_FD              = escorc(tmp_FD(3,:), tmp_FD(0,:))   ; Pearson’s correlation coefficient (r)

FD_GRACE_with_data= where( .not. ismissing(tmp_FD(0,:)), tmp_FD(1,:), default_value)
FD_GRACE_with_data@_FillValue = default_value
SD_FD             = stddev(FD_GRACE_with_data) ;stddev(tmp_FD(3,:))  ; SM
MBE_FD            = avg(tmp_FD(1,:) - tmp_FD(0,:))     ; Mean Bias Error (MBE)
tmp_all(3,:)      = tmp_FD(3,:)
delete(tmp_FD)

plots(0) = gsn_csm_xy (wks,x,tmp_all,res)  ; Create another filled XY plot.

; ________________ add order name ________________
txres                  = True
txres@txPerimOn        = False
txres@txFont           = 25
txres@txFontHeightF    = 0.025
amres                  = True
amres@amParallelPosF   = -0.487    ; This is the right edge of the plot.
amres@amOrthogonalPosF = -0.46    ; This is the bottom edge of the plot.
amres@amJust           = "TopLeft"

txid1   = gsn_create_text(wks, order_name(0), txres)
annoid1 = gsn_add_annotation(plots(0), txid1, amres)
; delete(amres@amJust)

; ________________ add r & MBE ________________
context = "~C~GRACE SD="+sprintf("%4.2f",SD_GRACE)  \
        + "~C~GW        SD="+sprintf("%4.2f",SD_GW) + ", r="+sprintf("%4.2f",r_GW) \
        + "~C~FD         SD="+sprintf("%4.2f",SD_FD) + ", r="+sprintf("%4.2f",r_FD)
delete(txres@txFont)
txres@txFontHeightF         = 0.02
txid2 = gsn_create_text(wks, context, txres)
amres@amParallelPosF        = -0.43    ; This is the right edge of the plot.
amres@amOrthogonalPosF      = -0.46     ; This is the bottom edge of the plot.
amres@amJust                = "TopLeft"
annoid2 = gsn_add_annotation(plots(0), txid2, amres)

; ________________ add shading ________________
gres             = True
gres@gsFillColor = "cornsilk" ;"rosybrown" ;"azure"

xmin = (/(2001-2000)*12+1,(2017-2000)*12+1/)
xmax = (/(2010-2000)*12,(2020-2000)*12/)
ymin = (/-100,-100/)
ymax = (/ 200, 200/)

xbars = (/ (/xmin(0), xmax(0), xmax(0), xmin(0), xmin(0)/), \
            (/xmin(1), xmax(1), xmax(1), xmin(1), xmin(1)/) /)

ybars = (/ (/ymin(0),ymin(0),ymax(0),ymax(0),ymin(0)/), \
           (/ymin(1),ymin(1),ymax(1),ymax(1),ymin(1)/)  /)
counter = 0
do shadow_num = 0,1
  poly(counter) = gsn_add_polygon(wks,plots(0),xbars(shadow_num,:),ybars(shadow_num,:),gres)
  counter = counter + 1
end do

gsn_panel(wks,(/plots/),(/1,1/),pres)

end
