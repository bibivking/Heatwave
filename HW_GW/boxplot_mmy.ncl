; ********************************************************
;copy from "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ********************************************************

undef("boxplot_mmy")
function boxplot_mmy (wks:graphic, x[*]:numeric, y[*][*]:numeric, boxOpts:logical, plotres:logical)

begin

  numbox    = dimsizes(y(:,0))
  boxWidth  = boxOpts@boxWidth
  boxcolors = boxOpts@boxColors
  xAxis     = plotres@tmXBValues

  plot = create "plot" logLinPlotClass wks
    ; "vpHeightF"    : 0.50
    ; "vpWidthF"     : 0.80
    ; "vpXF"         : 0.10
    ; "vpYF"         : 0.80
    "trYMinF"      : plotres@trYMinF
    "trYMaxF"      : plotres@trYMaxF
    "trXMinF"      : min(xAxis)
    "trXMaxF"      : max(xAxis)
    "pmTickMarkDisplayMode": "Always"
    "tmXBMode"             : "Explicit"
    "tmXBValues"           : xAxis
    ; "tmXBLabels"           : plotres@tmXBLabels
    "pmTitleDisplayMode"   : "Always"                   ; allow titles
    ; "tiMainOn"             : True
    ;"tiMainString"         : title
  end create

  ; attsetvalues(plot,res)

  polyres                    = True                          ;Set up defaults
  polyres@gsLineThicknessF   = 3.                         ; thickness of lines
  polyres@gsLineDashPattern  = 0

  ; Drawing boxplot
  dum = new((/numbox,9/),graphic)
  do gg = 0,numbox-1

    ff = xAxis(gg+1)
    polyres@gsLineColor = boxcolors(gg)
    ;
    ; ; Maximum horizontal solid line
    ; yy = (/y(gg,4),y(gg,4)/)
    ; xx = (/(ff-(boxWidth/4.)),(ff+(boxWidth/4.))/)
    ; polyres@gsLineDashPattern = 0
    ; if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
    ;      dum(gg,0) = gsn_add_polyline(wks,plot,xx,yy,polyres)
    ; end if
    ;
    ; ; 75th to Maximum vertical dash line
    ; yy = (/y(gg,3),y(gg,4)/)
    ; xx = (/ff,ff/)
    ; polyres@gsLineDashPattern = 1
    ; if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
    ;      dum(gg,1) = gsn_add_polyline(wks,plot,xx,yy,polyres)
    ; end if

    ; 75th horizontal solid line
    polyres@gsLineDashPattern = 0
    yy = (/y(gg,3),y(gg,3)/)
    xx = (/(ff-(boxWidth/2.)),(ff+(boxWidth/2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
         dum(gg,2) = gsn_add_polyline(wks,plot,xx,yy,polyres)
    end if

    ; 25th to 75th left vertical solid line
    yy = (/y(gg,1),y(gg,3)/)
    xx = (/(ff-(boxWidth/2.)),(ff-(boxWidth/2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
         dum(gg,3) = gsn_add_polyline(wks,plot,xx,yy,polyres)
    end if

    ; mean horizontal solid line
    yy = (/y(gg,2),y(gg,2)/)
    xx = (/(ff-(boxWidth/2.)),(ff+(boxWidth/2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
         dum(gg,4) = gsn_add_polyline(wks,plot,xx,yy,polyres)
    end if

    ; 25th to 75th right vertical solid line
    yy = (/y(gg,1),y(gg,3)/)
    xx = (/(ff+(boxWidth/2.)),(ff+(boxWidth/2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
         dum(gg,5) = gsn_add_polyline(wks,plot,xx,yy,polyres)
    end if

    ; 25th horizontal solid line
    yy = (/y(gg,1),y(gg,1)/)
    xx = (/(ff-(boxWidth/2.)),(ff+(boxWidth/2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
         dum(gg,6) = gsn_add_polyline(wks,plot,xx,yy,polyres)
    end if
    ;
    ; ; Minimum to 25th vertical dash line
    ; yy = (/y(gg,0),y(gg,1)/)
    ; xx = (/ff,ff/)
    ; polyres@gsLineDashPattern = 1
    ; if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
    ;      dum(gg,7) = gsn_add_polyline(wks,plot,xx,yy,polyres)
    ; end if
    ; polyres@gsLineDashPattern = 0
    ;
    ; ; Minimum horizontal solid line
    ; yy = (/y(gg,0),y(gg,0)/)
    ; xx = (/(ff-(boxWidth/4.)),(ff+(boxWidth/4.))/)
    ; if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
    ;      dum(gg,8) = gsn_add_polyline(wks,plot,xx,yy,polyres)
    ; end if

  end do

  dumname = unique_string("dum")
  plot@$dumname$ = dum

  ; res = True
  ; res@gsnDraw      = False
  ; res@gsnFrame     = False
  ; res@gsnMaximize  = False
  ;
  ; res@txFont                 = 25
  ; res@txFontHeightF          = 0.015
  ; res@tmBorderThicknessF     = 3.
  ;
  ; ; Y axis
  ; res@tmYLOn                 = True
  ; res@tmYLBorderOn           = True
  ; res@tmYLLabelsOn           = True
  ; res@tmYUseLeft             = True
  ; res@tmYLMajorThicknessF    = 3.
  ; res@tmYLMinorThicknessF    = 3.
  ; res@tmYLLabelFontThicknessF= 3.
  ; res@tmYLLabelFontHeightF   = 0.05
  ; res@tmYROn                 = True
  ; res@tmYRBorderOn           = False
  ; res@tmYRLabelsOn           = True
  ; res@tiYAxisFontThicknessF  = 3.
  ; res@tiYAxisFontHeightF     = 0.05
  ;
  ; ; X axis
  ; res@tmXBOn                 = True
  ; res@tmXBBorderOn           = True
  ; res@tmXBLabelsOn           = True
  ; res@tmXUseBottom           = True
  ; res@tmXBMajorThicknessF    = 3.
  ; res@tmXBMinorThicknessF    = 3.
  ; res@tmXBLabelFontThicknessF= 3.
  ; res@tiXAxisFontThicknessF  = 3.
  ; res@tmXBLabelFontHeightF   = 0.05
  ; res@tiXAxisFontHeightF     = 0.05
  ; res@tmXTOn                 = False
  ; res@tmXTBorderOn           = False
  ; res@tmXTLabelsOn           = False
  ; res@tmXBMode               = "Explicit"
  ; res@tmXBValues             = plotres@tmXBValues
  ; res@tmXBLabels             = plotres@tmXBLabels
  ;
  ; res@pmTitleDisplayMode     = "Always"
  ; res@pmLegendWidthF         = 0.2                     ; mresize legend width
  ; res@pmLegendHeightF        = 0.2                     ; and height
  ;
  ; ; res@gsnPaperOrientation = "landscape";"portrait";"landscape" ;"auto"
  ; maximize_output(wks,res)

  return(plot)
end
