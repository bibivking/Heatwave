;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./boxplot_mmy_old.ncl"
load "./boxplot_stats_mmy.ncl"
;*******************************************

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

; ================================== Options ===================================
year_s       = 2000
year_e       = 2019
operation    = "plot" ;"calc" ; "plot"
scale        = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
heatwave     = "hw-90th" ; "hw-90th" ; "all-day"
case_name    = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr", \
                 "CTL-30x3+30yr-DeepRoot/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr-DeepRoot"/)
message      = "CTL_vs_FD_vs_DR"
pic_name     = "HW_box-pdf_plot_"+heatwave+"_"+scale+"_"+message

var_name     = (/"Qle", "Qh", "Fwsoil", "deltaT", "SM64cm"/)

order_name   = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)

; ================================= Predefine ==================================
case_sum  = dimsizes(case_name)
year_sum = year_e-year_s+1

day_sum  = 0
do year = year_s,year_e
  if ( mod(year,4) .eq. 0) then
     day_sum = day_sum + 366
  else
     day_sum = day_sum + 365
  end if
end do

; _________ 4 time periods __________
          ;  all years,  2001-2009,  2017-2019,   2010-2011
          ;   2000-1-1,   2001-1-1,   2017-1-1,    2010-1-1
ts_s    = (/         0,        366,       6210,       3653/)
          ; 2019-12-31, 2009-12-31, 2019-12-31,  2011-12-31
ts_e    = (/ day_sum-1,       3652,       7304,       4382/)

total_box  = case_sum*2*4+3 ; case_sum x PFT all/tree x 4 time periods + space
bins       = 100
nVar       = case_sum*8

wks = gsn_open_wks("pdf","./plots/"+pic_name)

do var_num = 0, 4

  ; =============================== READ IN ======================================
  txt_name   = "./txt/"+var_name(var_num)+"_box_stats.txt"
  data_stats = asciiread(txt_name,(/total_box,5/),"float")
  txt_name   = "./txt/"+var_name(var_num)+"_pdf.txt"
  data_pdf   = asciiread(txt_name,(/nVar+1,bins/),"float")

  ; ================================ Plotting ====================================

  if var_name(var_num) .eq. "Qle" then
    ranges = (/0., 212./)
    bin_bdy= (/0., 250./)
    y_bdy  = (/0., 10./)
    y_title = "LH (W/m2)"
  else if var_name(var_num) .eq. "Qh" then
    ranges = (/0., 212./)
    bin_bdy= (/-50., 350./)
    y_bdy  = (/0., 2.9/)
    y_title = "SH (W/m2)"
  else if var_name(var_num) .eq. "Fwsoil" then
    ranges = (/0., 1.0/)
    bin_bdy= (/0., 1./)
    y_bdy  = (/0., 2.9/)
    y_title = "Fwsoil (-)"
  else if var_name(var_num) .eq. "deltaT" then
    ranges = (/0., 5.6/)
    bin_bdy= (/-3., 6./)
    y_bdy  = (/0., 5.4/)
    y_title = "Tleaf - Tair (C)"
  else if var_name(var_num) .eq. "SM64cm" then
    ranges = (/0.1, 0.312/)
    bin_bdy= (/0.05, 0.45/)
    y_bdy  = (/0., 5./)
    y_title = "SM (m3/m3)"
  end if
  end if
  end if
  end if
  end if

  x   = ispan(1,total_box,1)

  res            = True
  ; res@trYMinF    = ranges(0)
  ; res@trYMaxF    = ranges(1)
  ; res@tmXBValues = ispan(0,total_box+1,1)
  ; res@tmXBLabels = (/"","2000-19","","","", "T","","","", "",  \
  ;                       "2001-09","","","", "T","","","", "",  \
  ;                       "2017-19","","","", "T","","","", "",  \
  ;                       "2010-11","","","", "T","","","", ""   /)
  res@tmXBLabels = (/   "","","","","2000-19","","","","",  \
                        "","","","","2001-09","","","","",  \
                        "","","","","2010-11","","","","",  \
                        "","","","","2017-19","","",""      /)
  llres                   = True
  llres@gsLineThicknessF  = 3.                ; line thickness


  ; resources that control color and width of boxes
  opti           = True
  opti@boxWidth  = 0.9
  opti@boxColors = (/ \
                     "blue","red","green","orange","blue","red","green","orange","white", \
                     "blue","red","green","orange","blue","red","green","orange","white", \
                     "blue","red","green","orange","blue","red","green","orange","white", \
                     "blue","red","green","orange","blue","red","green","orange"\
                     /)

  ; plot = boxplot_mmy(wks, x, data_stats, opti, res)
  plot = boxplot_mmy(wks,x,data_stats,opti,res,llres, ranges(0), ranges(1), y_title)

  if var_num .eq. 0 then
    textres = True
    gres    = True

    lineres = True
    lineres@lgLineColors      = (/"white","blue","red","green","orange"/) ;,"white","blue","red","green","orange"/) ; line colors
    lineres@lgLineThicknesses = 2.5                        ; line thicknesses
    lineres@LineLengthPercent = 4.                         ; expressed as %, 0->100, length of line
    gres@YPosPercent          = 85.  ;expressed as %, 0->100, sets position of top border of legend (Default = 95.)

    gres@XPosPercent          = 30.  ;expressed as %, 0->100, sets position of left border of legend (Default = 5.)
    lineres@lgDashIndexes     = (/5,5,5,5,5/)
    textres@lgLabels          = (/"non-Tree", "CTL", "FD", "DR", "FD+DR"/)
    plot = simple_legend(wks,plot,gres,lineres,textres)
    gres@XPosPercent          = 50.  ;expressed as %, 0->100, sets position of left border of legend (Default = 5.)
    lineres@lgDashIndexes     = (/0,0,0,0,0/)
    textres@lgLabels          = (/"Tree", "CTL", "FD", "DR", "FD+DR"/)
    plot = simple_legend(wks,plot,gres,lineres,textres)
  end if
  ;
  ; draw(wks)                                    ; boxplot does not call these
  ; frame(wks)                                   ; for you

  res1                 = True
  res1@gsnDraw         = False  ; Don't draw plots
  res1@gsnFrame        = False

  res1@tmXBBorderOn    = False
  res1@tmXBOn          = False
  res1@tmXUseBottom    = False
  res1@tmXTLabelsOn    = True

  res1@tmYLOn          = True
  res1@tmYLBorderOn    = True
  res1@tmYRBorderOn    = False
  res1@tmYROn          = False
  res1@tmYUseLeft      = False
  res1@trYReverse      = True

  ; res1@tiXAxisFontThicknessF = 3.
  res1@tiXAxisFontHeightF    = 0.01
  res1@tmXTLabelFontHeightF  = 0.01
  ; res1@tmXTMajorThicknessF   = 3.
  ; res1@tmXTMinorThicknessF   = 3.
  ; res1@tiYAxisFontThicknessF = 3.
  res1@tiYAxisFontHeightF    = 0.01
  res1@tmYLLabelFontHeightF  = 0.01
  res1@tmYRLabelFontHeightF  = 0.01
  ; res1@tmYLMajorThicknessF   = 3.
  ; res1@tmYLMinorThicknessF   = 3.
  ; res1@tmYRMajorThicknessF   = 3.
  ; res1@tmYRMinorThicknessF   = 3.

  ; res1@gsnXYBarChart            = True              ; Create bar plot
  ; res1@gsnXYBarChartOutlineOnly = True

  res1@xyLineThicknesses      = (/1.5,1.5,1.5,1.5/)
  res1@xyLineOpacityF         = 1.
  res1@xyLineColors           = (/"blue","red","green","orange"/)
  res1@xyDashPatterns         = (/ 0,0,0,0  /)

  res1@pmLegendDisplayMode    = "Always"            ; turn on legend
  res1@pmLegendSide           = "Right"             ; Change location of
  res1@pmLegendParallelPosF   = .25                 ; move units right
  res1@pmLegendOrthogonalPosF = -0.4                ; move units down
  res1@pmLegendWidthF         = 0.125               ; Change width and
  res1@pmLegendHeightF        = 0.15                ; height of legend.
  res1@lgPerimOn              = False                ; turn off/on box around
  res1@lgLabelFontHeightF     = .015                ; label font height
  res1@tiYAxisString          = "PDF (%)"
  res1@tiXAxisSide            = "Top"

  res1@trYMinF        = y_bdy(0)
  res1@trYMaxF        = y_bdy(1)
  res1@trXMinF        = bin_bdy(0)
  res1@trXMaxF        = bin_bdy(1)

  res1@pmLegendDisplayMode   = "Never"            ; turn on legend

  ; 2000-2019
  res1@tiXAxisString   = y_title

  res1@vpWidthF        = 0.17
  res1@vpHeightF       = 0.16
  res1@vpXF            = 0.20
  res1@vpYF            = 0.865
  plot1 = gsn_csm_xy (wks, data_pdf(0,:), data_pdf(5:8,:), res1)
  delete(res1@tiYAxisString)

  ; 2001-2009
  res1@tmYLOn          = False
  res1@tmYLBorderOn    = False
  res1@vpXF            = 0.376
  res1@vpYF            = 0.865
  plot2 = gsn_csm_xy (wks, data_pdf(0,:), data_pdf(13:16,:), res1)

  ; 2010-2011
  res1@vpXF            = 0.552
  res1@vpYF            = 0.865
  plot3 = gsn_csm_xy (wks, data_pdf(0,:), data_pdf(21:24,:), res1)

  ; 2017-2019
  res1@pmLegendDisplayMode     = "Never"            ; turn on legend
  res1@tmYROn          = True
  res1@tmYRLabelsOn    = True
  res1@tmYRBorderOn    = True
  res1@vpXF            = 0.728
  res1@vpYF            = 0.865
  plot4 = gsn_csm_xy (wks, data_pdf(0,:), data_pdf(29:32,:), res1)

  pres                      =  True
  pres@gsLineThicknessF     =  2.0                      ;-- polyline thickness
  pres@gsLineColor          = "blue"                    ;-- polyline color

  txres                   = True
  txres@txPerimOn         = False
  txres@txFont            = 25
  txres@txFontHeightF     = 0.023

  amres                   = True
  amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
  amres@amOrthogonalPosF  = -0.57    ; This is the bottom edge of the plot.
  amres@amJust            = "TopLeft"

  txid1   = gsn_create_text(wks, order_name(var_num), txres)
  annoid1 = gsn_add_annotation(plot, txid1, amres)

  draw(wks)                                    ; boxplot does not call these
  frame(wks)                                   ; for you
end do ;var_num
end
