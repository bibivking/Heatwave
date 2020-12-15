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
case_name   = (/"CTL","FREE_DRAIN"/) ;,"DEEP_ROOT"/)
path        = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
pic_message = "CTL_vs_FREE_DRAIN"

year_s         = 2000
year_e         = 2019

year_sum       = year_e-year_s+1
total          = year_sum*12
case_sum       = dimsizes(case_name)
default_value  = -99999.

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)

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

; set PFT
if PFT .eq. "all" then
  iveg_s = 0
  iveg_e = 20
else if PFT .eq. "tree" then
  iveg_s = 1
  iveg_e = 4
end if
end if


; ================================= PLOT ======================================
pic = "./plots/time_serial_GRACE_JPL06-"+grace_version+"_vs_"+cable_version+\
      "_PFT-"+PFT+"_"+pic_message+"_"+scale+"_"+time_scale
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")

res                        = True              ; plot mods desired
res@gsnDraw                = False             ; don't draw yet
res@gsnFrame               = False             ; don't advance frame yet
res@gsnMaximize            = True
res@vpHeightF 	           = 0.4               ; change aspect ratio of plot
res@vpWidthF 	             = 0.8
res@txFontHeightF          = 0.015            ; change title font heights
res@tiYAxisFontThicknessF  = 3.
res@tmXBMajorThicknessF    = 3.
res@tmXBMinorThicknessF    = 3.
res@tmYLMajorThicknessF    = 3.
res@tmYLMinorThicknessF    = 3.
res@tmYLLabelFontThicknessF= 3.
res@tmBorderThicknessF     = 3.

res@txFont                 = "times-bold"
res@xyMonoLineColor        = False             ; want colored lines
res@tmXBMode               = "Explicit"
res@tiXAxisString          = "Year"
res@tiYAxisString          = "TWSA (km/m3)"
res@tfPolyDrawOrder        = "PreDraw"

; modify legend
; res@pmLegendDisplayMode    = "Always"
res@lgPerimOn              = False      ; No legend perimeter.
res@lgJustification        = "TopRight"
; res@pmLegendSide           = "Bottom"
res@pmLegendParallelPosF   = 0.5       ; Move to right
res@pmLegendOrthogonalPosF = 0.5     ; Move into plot
res@pmLegendWidthF         = 0.3       ; Decrease width
res@pmLegendHeightF        = 0.1       ; Decrease height
res@lgBoxMinorExtentF      = 0.2       ; Shorten the legend lines
res@lgItemOrder            = (/ 5,4,3,2,1,0 /)   ; Reorder the legends

LineColors       = (/"navy", "darkgreen","cornflowerblue","peru","darkolivegreen1","orchid"/) ; colors chosen
LineThicknesses	 = (/2.5,2.5,1.5,1.5,1.5,1.5/)      ; line thicknesses
DashPatterns	   = (/0.,0.,0.,0.,0.,0./)      ; make all lines solid
Opacities        = (/0.9,0.9,0.9,0.9,0.9,0.9/)
Labels           = (/"GRACE","TWSA", "AM","SM", "CInt", "SWE" /)

res@xyLineColors           = LineColors
res@xyLineThicknesses      = LineThicknesses
res@xyDashPatterns	       = DashPatterns
res@xyExplicitLegendLabels = Labels
res@xyLineOpacities        = Opacities
x                          = ispan(1,total,1)
res@tmXBValues             = ispan(1,total+1,12*2)
res@tmXBLabels             = tostring(ispan(year_s,year_e+1,2))
res@trYMaxF                = 200.
res@trYMinF                = -100.

gres             = True
gres@gsFillColor = "cornsilk" ;"rosybrown" ;"azure"

pres              = True
pres@wkWidth      = 2500
pres@wkHeight     = 2500
pres@gsnMaximize  = True

txres               = True
txres@txPerimOn     = False
; txres@txFont        = "times-bold"
txres@txFontHeightF = 0.018

amres               = True

plots = new(case_sum,graphic)
poly  = new(case_sum*2,graphic)

do case_num = 0,case_sum-1
  ; read ascii
  txt_name = case_name(case_num)+"_GRACE_"+year_s+"-"+year_e+".txt"
  tmp      = asciiread(txt_name,(/6,total/),"float")
  tmp@_FillValue = default_value
  r        = escorc(tmp(1,:), tmp(0,:))   ; Pearson’s correlation coefficient (r)
  MBE      = avg(tmp(1,:) - tmp(0,:))     ; Mean Bias Error (MBE)
  plots(case_num) = gsn_csm_xy (wks,x,tmp,res)  ; Create another filled XY plot.

  ; add order
  if case_num .eq. 0 then
    txid1 = gsn_create_text(wks, "(a)", txres)
    res@pmLegendDisplayMode    = "Always"
  else if case_num .eq. 1 then
    txid1 = gsn_create_text(wks, "(b)", txres)
    res@pmLegendDisplayMode    = "Never"
  else if case_num .eq. 2 then
    txid1 = gsn_create_text(wks, "(c)", txres)
    res@pmLegendDisplayMode    = "Never"
  end if
  end if
  end if
  amres@amParallelPosF        = -0.485    ; This is the right edge of the plot.
  amres@amOrthogonalPosF      = -0.46    ; This is the bottom edge of the plot.
  amres@amJust                = "TopLeft"
  annoid1 = gsn_add_annotation(plots(case_num), txid1, amres)
  delete(amres@amJust)

  ; add r & MBE
  context = "r        =  "+sprintf("%4.2f",r) + "~C~MBE = "+sprintf("%4.2f",MBE); ~C~ change line
  txid2 = gsn_create_text(wks, context, txres)
  amres@amParallelPosF        = -0.3    ; This is the right edge of the plot.
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
gsn_panel(wks,(/plots/),(/1,case_sum/),pres)

end
