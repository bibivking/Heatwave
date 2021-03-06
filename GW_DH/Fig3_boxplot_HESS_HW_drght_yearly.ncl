
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
load "./Fig3_boxplot_draw_mmy.ncl"
load "./Fig3_boxplot_stats_mmy.ncl"
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
case_name    = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr"/)
message      = "CTL_vs_FD"
pic_name     = "HW_boxplot_yearly_tree_"+heatwave+"_"+scale+"_"+message

var_name     = (/"deltaT", "EF", "Fwsoil", "TVeg", "GPP"/)
var_sum      = dimsizes(var_name)
order_name   = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)

; ================================= Ptomatoefine ==================================
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

; _________ time periods __________
        ; 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,\
        ; 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019
ts_s    = (/335, 700, 1065, 1430, 1796, 2161, 2526, 2891, 3257, 3622,\
           3987,4352, 4718, 5083, 5448, 5813, 6179, 6544, 6909/)

ts_e    = (/425, 790, 1155, 1521, 1886, 2251, 2616, 2982, 3347, 3712, \
          4077, 4443, 4808, 5173, 5538, 5904, 6269, 6634, 6999/)

total_box  = case_sum*2*dimsizes(ts_s)+dimsizes(ts_s)-1 ; case_sum x PFT all/tree x 4 time periods + space


; =============================== Pre-setting ==================================
wks = gsn_open_wks("pdf","./plots/Fig3_"+pic_name)

total_box1  = case_sum*dimsizes(ts_s)+dimsizes(ts_s)-1
x           = ispan(1,total_box1,1)

; ____________ setting for boxes ____________
res                = True

res@tmXBLabels = (/  \
                      "2001","","", "","","",     \
                      "2003","","", "","","",     \
                      "2005","","", "","","",     \
                      "2007","","", "","","",     \
                      "2009","","", "","","",     \
                      "2011","","", "","","",     \
                      "2013","","", "","","",    \
                      "2015","","", "","","",    \
                      "2017","","", "","","",    \
                      "2019","" /)

XBLabels = (/   "2001","2003","2005","2007","2009","2011",\
                "2013","2015","2017","2019" /)
xAxis_major    = (/1.5,7.5,13.5,19.5,25.5,31.5,37.5,43.5,49.5,55.5/)
llres                   = True
llres@gsLineThicknessF  = 3.                ; line thickness

; resources that control color and width of boxes
opti           = True
opti@boxWidth  = 0.9
opti@boxColors = (/ \
                   "dodgerblue2","tomato", "white","dodgerblue2","tomato", "white",\
                   "dodgerblue2","tomato", "white","dodgerblue2","tomato", "white",\
                   "dodgerblue2","tomato", "white","dodgerblue2","tomato", "white",\
                   "dodgerblue2","tomato", "white","dodgerblue2","tomato", "white",\
                   "dodgerblue2","tomato", "white","dodgerblue2","tomato", "white",\
                   "dodgerblue2","tomato", "white","dodgerblue2","tomato", "white",\
                   "dodgerblue2","tomato", "white","dodgerblue2","tomato", "white",\
                   "dodgerblue2","tomato", "white","dodgerblue2","tomato", "white",\
                   "dodgerblue2","tomato", "white","dodgerblue2","tomato", "white",\
                   "dodgerblue2","tomato" \
                   /)


; ____________ setting for diff lines ____________
res1 = res
res1@vpXF = 0.15
res1@vpYF = 0.70
res1@vpHeightF = 0.40 ;for 2 droughts ;0.50
res1@vpWidthF  = 0.70

; res1@tmXBBorderOn    = False
; res1@tmXBOn          = False
; res1@tmXUseBottom    = False
; res1@tmXTLabelsOn    = True
;
res1@tmYROn          = True
res1@tmYLOn          = False
; res1@tmYLBorderOn    = True
; res1@tmYRBorderOn    = False

res1@tmYUseLeft        = False
res1@tmYRLabelsOn      = True
; res1@tmYRMode          = "Explicit"
; res1@tmYRValues        = yrvalues
; res1@tmYRLabels        = y2labels
res1@tmXBOn             = False
res1@tmXTBorderOn       = True;False
res1@tmXTLabelsOn       = False;False
; res1@gsnPaperOrientation= "landscape"

res1@trXMinF            = 0
res1@trXMaxF            = 57
res1@tmXBLabelFontHeightF = 0.015
res1@tmYRLabelFontHeightF = 0.018
res1@xyLineThicknessF   = 3.0         ; make second line thicker
res1@xyDashPattern      = 0
res1@xyLineColor        = "gray40"
res1@tiYAxisSide        = "Right"
res1@tiYAxisFontHeightF = 0.023


; ____________ order text ____________
txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.023

amres                   = True
amres@amParallelPosF    = 0.42   ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.45;-0.88    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

; ____________ gsn_panel ____________
pres                            = True
pres@vpXF                       = 0.1
pres@vpYF                       = 0.70
pres@wkWidth                    = 2500
pres@wkHeight                   = 2500
pres@gsnMaximize                = False
pres@gsnPanelYWhiteSpacePercent = 3
; pres@gsnPanelXWhiteSpacePercent = 1

; _________ add shading ____________
gres             = True
gres@gsFillColor = "cornsilk" ;"rosybrown" ;"azure"

; ________ create graphic ___________
plots   = new(3,graphic)
plots1  = new(3,graphic)

plots_s  = new(3,graphic)
plots_s1 = new(3,graphic)

poly     = new(2,graphic)

; ____________ loop var ____________
do var_num = 0, var_sum -1

  if var_name(var_num) .eq. "deltaT" then
    ranges      = (/-0.5, 5/) ;(/-0.5, 5/)
    ranges_diff = (/-0.8, 0./)
    y_title = "~F5~D"+"~F25~T (~S~o~N~C)"
    y_title_R = "~F5~D"+"~F25~T~B~_GW~N~ - "+\
                "~F5~D"+"~F25~T~B~_FD~N~ (~S~o~N~C)"
  end if
  if var_name(var_num) .eq. "EF" then
    ranges      = (/0., 0.8/)
    ranges_diff = (/0., 0.2/)
    y_title = "~F25~EF (-)"
    y_title_R = "~F5~D"+"~F25~EF (-)"
  end if
  if var_name(var_num) .eq. "Fwsoil" then
    ranges      = (/0., 1.05/)
    ranges_diff = (/0., 0.4/)
    y_title = "~F5~b"+"~F25~ (-)"
    y_title_R = "~F5~Db"+"~F25~ (-)"
  end if
  if var_name(var_num) .eq. "GPP" then
    ranges      = (/0., 12/)
    ranges_diff = (/0., 1.5/)
    y_title = "~F25~GPP"+"~F25~ (gC d~S~-1~N~)"
    y_title_R = "~F5~D"+"~F25~GPP (gC d~S~-1~N~)"
  end if
  if var_name(var_num) .eq. "TVeg" then
    ranges      = (/0., 3.8/)
    ranges_diff = (/0., 1.1/)
    y_title = "~F25~Et (mm d~S~-1~N~)"
    y_title_R = "~F5~D"+"~F25~Et (mm d~S~-1~N~)"
  end if

  ; ____________ READ IN ____________
  txt_name       = "./txt/"+var_name(var_num)+"_CTL_FD_yearly_box_stats.txt"
  data_stats_tmp = asciiread(txt_name,(/total_box,5/),"float")
  data_stats     = new((/total_box1,5/),float, -9999.)
  do year = 0,  dimsizes(ts_s) - 1
    stempt1      = year*3
    stempt2      = year*5
    if year .lt. dimsizes(ts_s) - 1 then
      data_stats(stempt1:stempt1+2,:) = data_stats_tmp(stempt2+2:stempt2+4,:)
    else
      data_stats(stempt1:stempt1+1,:) = data_stats_tmp(stempt2+2:stempt2+3,:)
    end if
  end do
  delete(data_stats_tmp)

  ; ____________ Calc diff ____________
  median_diff = new(dimsizes(ts_s),float, -9999.)
  median_diff = data_stats(0::3,2) - data_stats(1::3,2) ; CTL - FD
  print(median_diff)
  median_x    = fspan(1.5,55.5,19)
  print(var_name)
  print(median_diff)

  ; ; ____________ Plotting shading ____________
  ; xmin = (/0 ,48/)
  ; xmax = (/27,57/)
  ; ymin = (/ranges(0),ranges(0)/)
  ; ymax = (/ranges(1),ranges(1)/)
  ;
  ; xbars = (/ (/xmin(0), xmax(0), xmax(0), xmin(0), xmin(0)/), \
  ; (/xmin(1), xmax(1), xmax(1), xmin(1), xmin(1)/) /)
  ;
  ; ybars = (/ (/ymin(0),ymin(0),ymax(0),ymax(0),ymin(0)/), \
  ; (/ymin(1),ymin(1),ymax(1),ymax(1),ymin(1)/)  /)
  ;
  ; if var_num .le. 2 then
  ;   do shadow_num = 0,1
  ;     poly(shadow_num) = gsn_add_polygon(wks,plots1(var_num),xbars(shadow_num,:),ybars(shadow_num,:),gres)
  ;   end do
  ; else
  ;   do shadow_num = 0,1
  ;     poly(shadow_num) = gsn_add_polygon(wks,plots_s1(var_num-3),xbars(shadow_num,:),ybars(shadow_num,:),gres)
  ;   end do
  ; end if

  ; ____________ Plotting box ____________
  if var_num .le. 2 then
    plots1(var_num)   = boxplot_mmy(wks,x,XBLabels,xAxis_major,data_stats,opti,res,llres, ranges(0), ranges(1), y_title)
  else
    plots_s1(var_num-3) = boxplot_mmy(wks,x,XBLabels,xAxis_major,data_stats,opti,res,llres, ranges(0), ranges(1), y_title)
  end if

  if var_num .eq. 0 .or. var_num .eq. 3 then
    textres = True
    gres    = True

    lineres = True
    lineres@lgLineColors      = (/"white","dodgerblue2","tomato","gray40"/);,"orange"/) ;,"white","dodgerblue2","tomato","green","orange"/) ; line colors
    lineres@lgLineThicknesses = 2.5                        ; line thicknesses
    lineres@LineLengthPercent = 4.                         ; expressed as %, 0->100, length of line

    gres@XPosPercent          = 25.  ;expressed as %, 0->100, sets position of left border of legend (Default = 5.)
    lineres@lgDashIndexes     = (/0,0,0/);,0/)
    textres@lgLabels          = (/" ", "GW", "FD", "GW - FD"/);, "FD+DR"/)
    if var_num .eq. 0 then
      gres@YPosPercent          = 25 ;30.  ;expressed as %, 0->100, sets position of top border of legend (Default = 95.)
      plots1(var_num)   = simple_legend(wks,plots1(var_num),gres,lineres,textres)
    else if var_num .eq. 3 then
      gres@YPosPercent          = 95 ;30.  ;expressed as %, 0->100, sets position of top border of legend (Default = 95.)
      plots_s1(var_num-3) = simple_legend(wks,plots_s1(var_num-3),gres,lineres,textres)
    end if
    end if
  end if

  if var_num .le. 2 then
    txid1   = gsn_create_text(wks, order_name(var_num), txres)
    annoid1 = gsn_add_annotation(plots1(var_num), txid1, amres)
  else
    txid1   = gsn_create_text(wks, order_name(var_num-3), txres)
    annoid1 = gsn_add_annotation(plots_s1(var_num-3), txid1, amres)
  end if

  res@tfDoNDCOverlay = True


  ; ____________ Plotting line ____________
  res1@tfDoNDCOverlay     = True
  res1@trYMinF            = ranges_diff(0)
  res1@trYMaxF            = ranges_diff(1)
  res1@tiYAxisString      = y_title_R
  res1@tmXTOn             = False
  if var_num .le. 2 then
    plots(var_num) =  gsn_csm_xy(wks,median_x,median_diff,res1)
    overlay(plots(var_num),plots1(var_num))
    draw(plots(var_num))
  else
    plots_s(var_num-3) =  gsn_csm_xy(wks,median_x,median_diff,res1)
    overlay(plots_s(var_num-3),plots_s1(var_num-3))
    draw(plots_s(var_num-3))
  end if

  frame(wks)

  ; plots(var_num) = plot
  ; delete(plot)
end do ;var_num

gsn_panel(wks,(/plots/),(/3,1/),pres)

gsn_panel(wks,(/plots_s/),(/3,1/),pres)

end
