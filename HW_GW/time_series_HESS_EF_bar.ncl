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

begin

; ____________________________ OPTIONS _______________________________
hydr_param = True
time_scale = "seasonal" ;"monthly" ; "daily"
scale      = "SE-AU"
year_s     = 2000
year_e     = 2019
year_sum   = year_e-year_s+1
time_tot   = year_sum*12

; ______________________________ Path _________________________________
path_mum       = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
case_name      = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr"/) ;,"CTL-30x3+30yr-DeepRoot/outputs"/)
pic_message    = "CTL_vs_FREE_DRAIN"
pic_name       = "HESS_EF_"+pic_message+"_"+time_scale+"_"+scale+"_"+year_s+"-"+year_e

order_name     = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)

; "tree" PFT
iveg_s         = 1
iveg_e         = 4

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

;______________________________ Data _________________________________
case_sum  = dimsizes(case_name)
path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = path_mum+case_name(case_num)
end do

; ___________________________ TREE Mask ______________________________
filename1   = path(0)+"/cable_out_1990_SE_Aus.nc"
f1          = addfile (filename1,"r")
PFT_mask    = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
default_val = f1->Qle@_FillValue

delete(filename1)
delete(f1)

; ___________________________ READ DATA ______________________________
Qle  = new((/case_sum,2,time_tot/),float,default_val)
Qh   = new((/case_sum,2,time_tot/),float,default_val)
Rnet = new((/case_sum,2,time_tot/),float,default_val)
Qle  = 0.
Qh   = 0.
Rnet = 0.

do case_num = 0, case_sum -1
  print("case_num = " +case_num)
  i = 0
  do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
    else
       dom = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
    end if

    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
    f1 = addfile (filename1,"r")
    do month_num = 0,11
      Qle(case_num, 0, i)  = avg(f1->Qle(dom(month_num):dom(month_num+1)-1,{A:B},{C:D}))
      Qh(case_num, 0, i)   = avg(f1->Qh(dom(month_num):dom(month_num+1)-1,{A:B},{C:D}))
      Rnet(case_num, 0, i) = avg(f1->Rnet(dom(month_num):dom(month_num+1)-1,{A:B},{C:D}))
      Qle(case_num, 1, i)  = avg(where(PFT_mask .eq. 1, \
                                       dim_avg_n_Wrap(f1->Qle(dom(month_num):dom(month_num+1)-1,{A:B},{C:D}), 0),\
                                       default_val))
      Qh(case_num, 1, i)   = avg(where(PFT_mask .eq. 1, \
                                       dim_avg_n_Wrap(f1->Qh(dom(month_num):dom(month_num+1)-1,{A:B},{C:D}), 0),\
                                       default_val))
      Rnet(case_num, 1, i) = avg(where(PFT_mask .eq. 1, \
                                       dim_avg_n_Wrap(f1->Rnet(dom(month_num):dom(month_num+1)-1,{A:B},{C:D}), 0),\
                                       default_val))
      i = i + 1
    end do ; month_num
  end do ; year
end do ; case_num

print("prepare to plot")

if time_scale .eq. "seasonal" then
  Qle!0 = "case"
  Qle!1 = "PFT"
  Qle!2 = "time"
  Qh!0  = "case"
  Qh!1  = "PFT"
  Qh!2  = "time"
  Qle_season  = new((/case_sum,2,year_sum*4/),float, default_val)
  Qh_season   = new((/case_sum,2,year_sum*4/),float, default_val)
  do case_num = 0, case_sum-1
    do pft = 0, 1
      print(dimsizes(month_to_seasonN(Qle(case_num,pft,:),(/"DJF","MAM","JJA","SON"/))))
      tmp_Qle = month_to_seasonN(Qle(case_num,pft,:),(/"DJF","MAM","JJA","SON"/))
      tmp_Qle!0 = "season"
      tmp_Qle!1 = "year"
      Qle_season(case_num,pft,:) = ndtooned(tmp_Qle(year|:,season|:))
      delete(tmp_Qle)
      tmp_Qh = month_to_seasonN(Qh(case_num,pft,:),(/"DJF","MAM","JJA","SON"/))
      tmp_Qh!0 = "season"
      tmp_Qh!1 = "year"
      Qh_season(case_num,pft,:)  = ndtooned(tmp_Qh(year|:,season|:))
      delete(tmp_Qh)
    end do
  end do
  printVarSummary(Qle_season)
end if
;_____________________________ PLOT ___________________________________
pic = "./plots/"+pic_name
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow") ;"BlueWhiteOrangeRed") ;"BlueYellowRed")

res                         = True
res@gsnDraw                 = False  ; Don't draw plots
res@gsnFrame                = False  ; ����ҳ
res@gsnMaximize             = True

res@vpHeightF               = 0.4
res@vpWidthF                = 0.8
; res@vpXF                    = 0.03
; res@vpYF                    = 0.7
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

; res@xyLineThicknesses       = (/ 1.5, 1.5, 1.5, 1.5/)
; res@xyDashPatterns          = (/   0,   0,   0,   0/)
; res@xyLineColors            = (/"blue","green","red","orange"/)
; res@xyExplicitLegendLabels  = (/"CTL_Qle/Rnet","FD_Qle/Rnet","CTL_Qh/Rnet","FD_Qh/Rnet"/)

; res@lbLabelBarOn            = True
; res@lbBoxLinesOn            = False                       ;�ر�lbar box ����
; res@lbTitleFont             = 25
; res@lbLabelFont             = 25
; res@lbTitleFontHeightF      = 0.013
; res@lbLabelFontHeightF      = 0.013

res@lgPerimOn              = False      ; No legend perimeter.
res@lgJustification        = "TopRight"
; res@pmLegendSide           = "Bottom"
res@pmLegendParallelPosF   = 0.5       ; Move to right
res@pmLegendOrthogonalPosF = 0.5     ; Move into plot
res@pmLegendWidthF         = 0.3       ; Decrease width
res@pmLegendHeightF        = 0.1       ; Decrease height
res@lgBoxMinorExtentF      = 0.2       ; Shorten the legend lines
res@lgItemOrder            = (/ 3,2,1,0 /)   ; Reorder the legends

res@tmXBLabelFont           = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont           = 25
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF    = 0.015
res@tmYLLabelFontHeightF    = 0.015

res@tmYLMode                = "Explicit"
res@tmXBMode                = "Explicit"
res@txString                = ""
res@gsnLeftString           = ""
res@gsnRightString          = ""
; res@tmYLValues          = fspan(-2.0,2.0,21)
; res@tmYLLabels          = fspan(-2.0,2.0,21)
; res@tiYAxisString       = "Depth (cm)"


;**************** ͼ�������� *****************
pres                    = True                                      ; ���� panel resources.
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.2
pres@lbLabelFontHeightF = 0.2

txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.03

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"


; ________________________________ Plot Qh Qle _________________________________
plots              = new(4,graphic)
res@gsnYRefLine    = 0.              ; reference line
res@gsnXYBarChart  = True            ; create bar chart
res@trYMinF        = -15.
res@trYMaxF        = 15.

res@gsnAboveYRefLineColor = "red"           ; above ref line fill red
res@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
res@pmLegendDisplayMode   = "Never"

if time_scale .eq. "seasonal" then
  time_tot        = year_sum*4
  res@tmXBValues  = ispan(1,time_tot,4*2)
  x               = ispan(1,time_tot,1)
  res@tmXBLabels  = ispan(year_s,year_e,2)
  res@tiXAxisString = "Year"
  plots(0)          = gsn_csm_xy(wks, x, Qle_season(1,0,:)-Qle_season(0,0,:),res)
  plots(1)          = gsn_csm_xy(wks, x, Qle_season(1,1,:)-Qle_season(0,1,:),res)
  plots(2)          = gsn_csm_xy(wks, x, Qh_season(1,0,:)-Qh_season(0,0,:),res)
  plots(3)          = gsn_csm_xy(wks, x, Qh_season(1,1,:)-Qh_season(0,1,:),res)
else if time_scale .eq. "monthly" then
  time_tot        = year_sum*12
  res@tmXBValues  = ispan(1,time_tot,12*2)
  x               = ispan(1,time_tot,1)
  res@tmXBLabels  = ispan(year_s,year_e,2)
  res@tiXAxisString = "Year"
  plots(0)          = gsn_csm_xy(wks, x, Qle(1,0,:)-Qle(0,0,:),res)
  plots(1)          = gsn_csm_xy(wks, x, Qle(1,1,:)-Qle(0,1,:),res)
  plots(2)          = gsn_csm_xy(wks, x, Qh(1,0,:)-Qh(0,0,:),res)
  plots(3)          = gsn_csm_xy(wks, x, Qh(1,1,:)-Qh(0,1,:),res)
end if
end if

do counter = 0,3
  txid1   = gsn_create_text(wks, order_name(counter), txres)
  annoid1 = gsn_add_annotation(plots(counter), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/2,2/),pres)

end
