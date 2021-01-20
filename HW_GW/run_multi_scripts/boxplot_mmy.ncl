; ********************************************************
;copy from "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ********************************************************

; Adam Phillips
;
undef("boxplot_mmy")
function boxplot_mmy (wks:graphic,x[*]:numeric,y[*][*]:numeric,boxOpts:logical,plotres:logical,lineres:logical,min_value:numeric,max_value:numeric,title:string)
;
;    This function creates a blank linLog plot object, on which box plots are created by extensive
;    use of gsn_add_polyline. The user can draw as many or as few box plots as they wish.
;    The plot is not drawn, and the frame is not advanced. May be used with gsn_panel. Missing data
; may be present in the input Y array, but not the input X array.
;    The function options are as follows:
;
;    wks     ----- The workstation
;    x[*]    ----- A one-dimensional array containing the X-axis values of where the box plot(s)
;                  will be drawn.
;    y[*][*] ----- A two-dimensional array, where the rightmost dimension contains the box plot
;               reference pt. data. y(n,0)=bottom_value, y(n,1)=bottom_value_of_box,
;               y(n,2)=mid-value_of_box,y(n,3)=top_value_of_box,y(n,4)=top_value
;    boxOpts ----- 2 options attached as attributes may be attached here.
;               boxWidth ---- Scalar or array that contains the widths of the boxes.
;               boxColors ---- Scalar or array that contains the colors that the boxes will be drawn in.
;                           Ex. boxOpts@boxColors =  (/"green","blue"/)
;                        If the number of colors specified does not equal the number of
;                        boxes drawn, only the first color specified will be used.
;    plotres ----- An optional xy-plot resource list. Will override defaults where applicable.
;    lineres ----- An optional resource list pertaining to the lines used to draw the boxes/lines.
;
;    Example call: plot3 = boxplot(wks,ntime,newyval,opti,res,False)
;               draw(plot3)
;               frame(wks)
;
begin
     dimquery = dimsizes(y) ; data dims
     numbox = dimquery(0)   ; box numbers
     boxWidths = new((/numbox/),float) ; set box widths
     boxWidths = boxOpts@boxWidth

     ; if (numbox.ne.dimsizes(x)) then
     ;    print("boxplot: Fatal: X must be one-dimensional and both X and Y must have the same leftmost dimension")
     ;    exit
     ; end if
     ; if (any(ismissing(x))) then
     ;    print("boxplot: Fatal: X array cannot contain missing data, exiting")
     ;    exit
     ; end if

     xAxis = new(numbox+2,typeof(x)) ; Developing x-axis
     xAxis(1:numbox) = x
     dx = x(1)-x(0) ; distance between two boxes
     xAxis(0) = x(0)-dx ; left boundary of the box
     xAxis(numbox+1) = x(numbox-1)+dx ; right boundary of the box

     labarr = new(numbox+2,"string")         ;Prepare actual X-axis labels...
     labarr(0) = ""
     labarr(numbox+1) = ""
     labarr(1:numbox) = xAxis(1:numbox)

; Whether to maximize plot in frame.
     ; maximize = get_res_value(plotres,"gsnMaximize",False)
;
;      if (plotres) then
; ;       print("Plot resources detected, accepting")
      ;   fsatts = getvaratts(plotres)
      ;   do ty = 0,dimsizes(fsatts)-1
      ;      if (fsatts(ty).eq."tmXBLabels") then        ;Special section to test for XBLabels
      ;         if (dimsizes(plotres@tmXBLabels).ne.numbox) then
      ;            print("boxplot: Fatal:Number of XB Labels does not match number of boxes, exiting")
      ;            exit
      ;         else
      ;            labarr(1:numbox) = plotres@$fsatts(ty)$
      ;         end if
      ;         xblab = plotres@tmXBLabels
      ;         delete(plotres@tmXBLabels)      ;Delete so tmXBLabels is not used when all
      ;      end if                       ;atts are assigned below...
      ;   end do
      ;   delete(fsatts)
      ; end if

      res = True
      res@gsnDraw= False
      res@gsnFrame= False
      res@vpHeightF= 0.50
      res@vpWidthF= 0.80
      res@vpXF= 0.10
      res@vpYF= 0.80

      res@tmXTBorderOn= False
      res@tmXTOn= False
      res@tmXUseBottom= True
      res@tmXTLabelsOn= False

      res@tmYRBorderOn= False
      res@tmYROn= False
      res@tmYUseLeft= False

      res@tiYAxisFontThicknessF= 3.
      res@tmXBMajorThicknessF= 3.
      res@tmXBMinorThicknessF= 3.
      res@tmYLMajorThicknessF= 3.
      res@tmYLMinorThicknessF= 3.
      res@tmYLLabelFontThicknessF= 3.
      res@tmBorderThicknessF= 3.
      res@txFont= 25
      res@txFontHeightF= 0.015
      res@gsnMaximize= True                   ; maximize plot in frame
      res@pmLegendWidthF= 0.2                     ; mresize legend width
      res@pmLegendHeightF= 0.2                     ; and height
      res@tiXAxisFontHeightF   = 0.05
      res@tiYAxisFontHeightF   = 0.05
      res@tmXBLabelFontHeightF = 0.05
      res@tmYLLabelFontHeightF = 0.05

      res@trYMinF=               min_value ;-2 ! MMY
      res@trYMaxF=               max_value ;+2 ! MMY
      res@trXMinF=               min(xAxis)
      res@trXMaxF=               max(xAxis)
      res@pmTickMarkDisplayMode= "Always"
      res@tmXBMode=              "Explicit"
      res@tmXBValues=            "xAxis"
      res@tmXBLabels=            "labarr"
      res@pmTitleDisplayMode=    "Always"
      res@tiMainString=          title

      polyres = True                          ;Set up defaults
      polyres@gsLineColor       = "black"                     ; color of lines
      polyres@gsLineThicknessF  = 1.5                          ; thickness of lines
      polyres@gsLineDashPattern = 0
      if (lineres) then
        fsatts = getvaratts(lineres)
        do ty = 0,dimsizes(fsatts)-1
           polyres@$fsatts(ty)$ = lineres@$fsatts(ty)$
        end do
      end if

      boxcolor = boxOpts@boxColors

      dum = new((/numbox,9/),graphic)
      do gg = 0,numbox-1

        ff = xAxis(gg+1)
        polyres@gsLineColor = boxcolor(gg)

        ; Maximum horizontal solid line
        ; yy = (/y(gg,4),y(gg,4)/)
        ; xx = (/(ff-(boxWidths(gg)/8.)),(ff+(boxWidths(gg)/8.))/)
        ; polyres@gsLineDashPattern = 0
        ; if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
        ;      dum(gg,0) = gsn_add_polyline(wks,plot,xx,yy,polyres)
        ; end if

        ; 75th to Maximum vertical dash line
        ; yy = (/y(gg,3),y(gg,4)/)
        ; xx = (/ff,ff/)
        ; polyres@gsLineDashPattern = 1
        ; if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
        ;      dum(gg,1) = gsn_add_polyline(wks,plot,xx,yy,polyres)
        ; end if

        ; 75th horizontal solid line
        polyres@gsLineDashPattern = 0
        yy = (/y(gg,3),y(gg,3)/)
        xx = (/(ff-(boxWidths(gg)/2.)),(ff+(boxWidths(gg)/2.))/)
        if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(gg,2) = gsn_add_polyline(wks,plot,xx,yy,polyres)
        end if

        ; 25th to 75th left vertical solid line
        yy = (/y(gg,1),y(gg,3)/)
        xx = (/(ff-(boxWidths(gg)/2.)),(ff-(boxWidths(gg)/2.))/)
        if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(gg,3) = gsn_add_polyline(wks,plot,xx,yy,polyres)
        end if

        ; mean horizontal solid line
        yy = (/y(gg,2),y(gg,2)/)
        xx = (/(ff-(boxWidths(gg)/2.)),(ff+(boxWidths(gg)/2.))/)
        if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(gg,4) = gsn_add_polyline(wks,plot,xx,yy,polyres)
        end if

        ; 25th to 75th right vertical solid line
        yy = (/y(gg,1),y(gg,3)/)
        xx = (/(ff+(boxWidths(gg)/2.)),(ff+(boxWidths(gg)/2.))/)
        if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(gg,5) = gsn_add_polyline(wks,plot,xx,yy,polyres)
        end if

        ; 25th horizontal solid line
        yy = (/y(gg,1),y(gg,1)/)
        xx = (/(ff-(boxWidths(gg)/2.)),(ff+(boxWidths(gg)/2.))/)
        if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(gg,6) = gsn_add_polyline(wks,plot,xx,yy,polyres)
        end if

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
        ; xx = (/(ff-(boxWidths(gg)/8.)),(ff+(boxWidths(gg)/8.))/)
        ; if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
        ;      dum(gg,8) = gsn_add_polyline(wks,plot,xx,yy,polyres)
        ; end if

      end do
      dumname = unique_string("dum")
      plot@$dumname$ = dum
      if (isvar("xblab").eq."True") then
        plotres@tmXBLabels = xblab      ;reassign XBLabels for possible future plots
      end if
      if (maximize) then
        mres = True
        mres@gsnDraw  = False
        mres@gsnFrame = False
        mres@vpHeightF= 0.50
        mres@vpWidthF = 0.80
        mres@vpXF     = 0.10
        mres@vpYF     = 0.80


        mres@tmXTBorderOn    = False
        mres@tmXTOn          = False
        mres@tmXUseBottom    = True
        mres@tmXTLabelsOn    = False

        mres@tmYRBorderOn    = False
        mres@tmYROn          = False
        mres@tmYUseLeft      = False

        mres@tiYAxisFontThicknessF = 3.
        mres@tmXBMajorThicknessF   = 3.
        mres@tmXBMinorThicknessF   = 3.
        mres@tmYLMajorThicknessF   = 3.
        mres@tmYLMinorThicknessF   = 3.
        mres@tmYLLabelFontThicknessF= 3.
        mres@tmBorderThicknessF     = 3.
        mres@txFont                 = 25
        mres@txFontHeightF        = 0.015
        mres@gsnMaximize          = True                   ; maximize plot in frame
        mres@pmLegendWidthF       = 0.2                     ; mresize legend width
        mres@pmLegendHeightF      = 0.2                     ; and height
        mres@tiXAxisFontHeightF   = 0.05
        mres@tiYAxisFontHeightF   = 0.05
        mres@tmXBLabelFontHeightF = 0.05
        mres@tmYLLabelFontHeightF = 0.05
        ;
        ; mres@gsnPaperOrientation = "landscape";"portrait";"landscape" ;"auto"
        ; _________________________
        maximize_output(wks,mres)
      end if
     return(plot)
end
