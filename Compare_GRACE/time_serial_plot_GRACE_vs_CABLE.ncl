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
case_name   = (/"CTL","FD","DR","FD+DR"/) ;,"DEEP_ROOT"/)
path        = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
pic_message = "CTL-FD-DR-FD+DR"

year_s         = 2000
year_e         = 2019

year_sum       = year_e-year_s+1
total          = year_sum*12
case_sum       = dimsizes(case_name)
default_value  = -99999.

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)
order_name     = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)

; ================================ Set Regrid ==================================
PFT           = "tree" ; "tree" ;"all"
scale         = "SE-AU"
time_scale    = "monthly" ; "annual"; "monthly"
cable_version = "CABLE-HESS" ;"LIS-CABLE" ; "CABLE-HESS"
grace_version = "Original-fill-nan"
                ;"Original-fill-nan": using Nan fill missing month
                ;"Gap-fill" : Sanaa's gap fill
                ;"Original" : if this option, should alwasy set GRACE_start_year = 2003
                ;             GRACE_end_year <= 2010, since too many gaps after 2010


; ================================= PLOT ======================================
pic = "./plots/time_serial_GRACE_JPL06-"+grace_version+"_vs_"+cable_version+\
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

res@txFont                 = "times-bold"
res@xyMonoLineColor        = False             ; want colored lines
res@tmXBMode               = "Explicit"
; res@tiXAxisString          = "Year"
res@tfPolyDrawOrder        = "PreDraw"

res@lgPerimOn              = False      ; No legend perimeter.
; res@lgJustification        = "TopRight"
; res@pmLegendSide           = "Bottom"
res@pmLegendParallelPosF   = 0.7       ; Move to right
res@pmLegendOrthogonalPosF = -0.6     ; Move into plot
res@pmLegendWidthF         = 0.15       ; Decrease width
res@pmLegendHeightF        = 0.15       ; Decrease height
res@lgBoxMinorExtentF      = 0.3       ; Shorten the legend lines
; res@lgItemOrder            = (/ 2,1,0 /)   ; Reorder the legends

res@lgItemOrder            = (/ 5,4,3,2,1,0 /)   ; Reorder the legends
LineColors       = (/ "black","darkgreen","cornflowerblue","peru","darkolivegreen1","orchid"/) ; "navy",colors chosen
LineThicknesses	 = (/2.5,2.5,1.5,1.5,1.5,1.5/)      ; line thicknesses
DashPatterns	   = (/0.,0.,0.,0.,0.,0./)      ; make all lines solid
Opacities        = (/0.9,0.9,0.9,0.9,0.9,0.9/)
Labels           = (/"GRACE","CABLE", "Aquifer","SM", "CInt", "SWE" /)

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

gres             = True
gres@gsFillColor = "cornsilk" ;"rosybrown" ;"azure"

pres              = True
pres@wkWidth      = 2500
pres@wkHeight     = 2500
pres@gsnMaximize  = False
pres@gsnPanelYWhiteSpacePercent = 0.1
pres@gsnPanelXWhiteSpacePercent = 0.1

txres               = True
txres@txPerimOn     = False
; txres@txFont        = "times-bold"
txres@txFontHeightF = 0.025

amres               = True

plots = new(case_sum,graphic)
poly  = new(case_sum*2,graphic)

do case_num = 0,case_sum-1
  ; add legend

  if mod(case_num,2) .eq. 0 then
    res@tiYAxisString          = "TWSA (mm)"
  else
    delete(res@tiYAxisString)
  end if

  ; read ascii
  txt_name = case_name(case_num)+"_GRACE_PFT-"+PFT+"_"+year_s+"-"+year_e+".txt"
  tmp      = asciiread(txt_name,(/6,total/),"float")
  tmp@_FillValue = default_value
  r        = escorc(tmp(1,:), tmp(0,:))   ; Pearson’s correlation coefficient (r)
  MBE      = avg(tmp(1,:) - tmp(0,:))     ; Mean Bias Error (MBE)
  SD_GRACE = stddev(tmp(0,:))
  ; tmp2     = where(.not. ismissing(tmp(0,:)), tmp(1,:), default_value)
  ; tmp2@_FillValue = default_value
  SD_CABLE = stddev(tmp(1,:))
  ; delete(tmp2)

  r_SM     = escorc(tmp(3,:), tmp(0,:))   ; Pearson’s correlation coefficient (r)
  SD_SM    = stddev(tmp(3,:))



  if case_num .eq. 0 then
    res@pmLegendDisplayMode    = "Always"
    context = "r                   =  "+sprintf("%4.2f",r) +" ("+sprintf("%4.2f",r_SM)+")"\
            + "~C~GRACE SD = "+sprintf("%4.2f",SD_GRACE)       \
            + "~C~CABLE  SD = "+sprintf("%4.2f",SD_CABLE) +" ("+ sprintf("%4.2f",SD_SM)+")"
            ; ~C~ change line
  else
    res@pmLegendDisplayMode    = "Never"
    context = "r                   =  "+sprintf("%4.2f",r) +" ("+sprintf("%4.2f",r_SM)+")"\
            + "~C~GRACE SD = "+sprintf("%4.2f",SD_GRACE)       \
            + "~C~CABLE  SD = "+sprintf("%4.2f",SD_CABLE) +" ("+ sprintf("%4.2f",SD_SM)+")"
    ; context = "r                   =  "+sprintf("%4.2f",r) \
    ;         + "~C~GRACE SD = "+sprintf("%4.2f",SD_GRACE)   \
    ;         + "~C~CABLE  SD = "+sprintf("%4.2f",SD_CABLE)
            ; ~C~ change line
  end if

  ; plots(case_num) = gsn_csm_xy (wks,x,tmp(0:2,:),res)  ; Create another filled XY plot.
  plots(case_num) = gsn_csm_xy (wks,x,tmp(:,:),res)  ; Create another filled XY plot.
  amres@amParallelPosF        = -0.485    ; This is the right edge of the plot.
  amres@amOrthogonalPosF      = -0.46    ; This is the bottom edge of the plot.
  amres@amJust                = "TopLeft"
  txid1   = gsn_create_text(wks, order_name(case_num), txres)
  annoid1 = gsn_add_annotation(plots(case_num), txid1, amres)
  delete(amres@amJust)

  ; add r & MBE
  txid2 = gsn_create_text(wks, context, txres)
  amres@amParallelPosF        = -0.4    ; This is the right edge of the plot.
  amres@amOrthogonalPosF      = -0.46     ; This is the bottom edge of the plot.
  amres@amJust                = "TopLeft"
  annoid2 = gsn_add_annotation(plots(case_num), txid2, amres)

  delete(txt_name)
  delete(tmp)
  delete(r)
  delete(MBE)
end do ; case_num

xmin = (/(2001-2000)*12+1,(2017-2000)*12+1/)
xmax = (/(2010-2000)*12,(2020-2000)*12/)
ymin = (/-100,-100/)
ymax = (/ 200, 200/)

xbars = (/ (/xmin(0), xmax(0), xmax(0), xmin(0), xmin(0)/), \
            (/xmin(1), xmax(1), xmax(1), xmin(1), xmin(1)/) /)

ybars = (/ (/ymin(0),ymin(0),ymax(0),ymax(0),ymin(0)/), \
           (/ymin(1),ymin(1),ymax(1),ymax(1),ymin(1)/)  /)
counter = 0
do case_num = 0, case_sum-1
  do shadow_num = 0,1
    poly(counter) = gsn_add_polygon(wks,plots(case_num),xbars(shadow_num,:),ybars(shadow_num,:),gres)
  counter = counter + 1
  end do
end do
gsn_panel(wks,(/plots/),(/2,2/),pres)

end
