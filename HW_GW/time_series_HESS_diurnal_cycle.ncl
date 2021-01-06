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


; ================================= Plotting ===================================
res1                    = True
res1@cnFillMode         = "RasterFill"            ; Raster Mode
res1@cnFillOn           = True                    ;��ɫ
res1@tmBorderThicknessF = 3.0

res1@gsnDraw            = True  ; Don't draw plots
res1@gsnFrame           = True  ; ����ҳ
res1@lbLabelBarOn       = True

;************** ����labelbar ***************
res1@lbBoxLinesOn       = True                       ;�ر�lbar box ����
res1@lbTitleFont        = 25
res1@lbLabelFont        = 25
res1@lbTitleFontHeightF = 0.013
res1@lbLabelFontHeightF = 0.013
res1@txString           = ""
res1@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
res1@tmYLLabelFont      = 25

;*************** ���õ�ֵ�� *****************
res1@cnLinesOn          = False                       ; �رյ�ֵ������
res1@cnLineColor        = "black"
res1@cnLineThicknessF   = 1.5
res1@cnLineLabelsOn     = False
res1@gsnMaximize        = True
res1@cnExplicitLabelBarLabelsOn = True   ;?
res1@cnLevelSelectionMode = "ExplicitLevels"

;************ ����������ֵ��Χ **************
res1@tmXBLabelFontThicknessF = 0.015
res1@tmYLLabelFontThicknessF = 0.015
res1@tmXBLabelFontHeightF = 0.015
res1@tmYLLabelFontHeightF = 0.015
res1@tmYLMode  = "Explicit"
res1@tmXBMode  = "Explicit"

res1@gsnLeftString       = ""
res1@gsnRightString      = ""
wks1 = gsn_open_wks("pdf","check")
gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed");"WhiteBlueGreenYellowRed");"GMT_haxby");"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")


; ____________________________ OPTIONS _______________________________
PFT        = "tree" ; "tree"
scale      = "SE-AU"
year_s     = 2019
year_e     = 2019

HW_day_s   = (/ 15, 22/)
HW_day_e   = (/ 19, 27/)
event_sum  = dimsizes(HW_day_s)
ts_tot     = 9
day_tot    = 0
do event_num = 0, event_sum-1
  day_tot = day_tot + (HW_day_e(event_num)-HW_day_s(event_num)+1)
end do

message    = "subday-24h"

pic_name   = "CABLE-HESS_"+message+"_PFT-"+PFT+"_"+scale+"_"+year_s

path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"

case_name = (/"CTL-30x3+30yr","FREE_DRAIN",\
              "CTL-30x3+30yr_2019-DeepRoot",\
              "CTL-30x3+30yr_2019-WetSoil"/)

case_brf  = (/"non-Tree CTL","      FD","      CTL+DR","     CTL+WetS", \
              "Tree CTL","   FD","   CTL+DR","   CTL+WetS"/)

order_name  = (/"(a)","(b)","(c)","(d)","(e)","(f)",\
                "(g)","(h)","(i)","(j)","(k)","(l)"/)

; 0, 1, 2, 3 , 4 , 5,  6,  7
; 3, 6, 9, 12, 15, 18, 21, 24

iveg_tree_s = 1
iveg_tree_e = 4
iveg_nontree_s = 5
iveg_nontree_e = 17

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
  path(case_num) = path_mum+case_name(case_num)+"/outputs-subday"
end do

; ___________________________ TREE Mask ______________________________
filename1 = path(0)+"/cable_out_2019_SE_Aus.nc"
f1 = addfile (filename1,"r")
PFT_mask = where((f1->iveg({A:B},{C:D}) .ge. iveg_tree_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_tree_e), 1, f1->iveg({A:B},{C:D}))
PFT_mask = where((PFT_mask .ge. iveg_nontree_s) .and. (PFT_mask .le. iveg_nontree_e), 2, PFT_mask)
default_val = f1->Qle@_FillValue

PFT_mask_multi = new((/ts_tot,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
do z = 0, ts_tot-1
  PFT_mask_multi(z,:,:) = PFT_mask
end do
plot = gsn_csm_contour(wks1,PFT_mask,res1)
delete(filename1)
delete(f1)

Tair     = new((/case_sum*2,day_tot,ts_tot/),float,default_val)
Qle      = new((/case_sum*2,day_tot,ts_tot/),float,default_val)
Qh       = new((/case_sum*2,day_tot,ts_tot/),float,default_val)
deltaT   = new((/case_sum*2,day_tot,ts_tot/),float,default_val)
Fwsoil   = new((/case_sum*2,day_tot,ts_tot/),float,default_val)
Rnet     = new((/case_sum*2,day_tot,ts_tot/),float,default_val)
Qg       = new((/case_sum*2,day_tot,ts_tot/),float,default_val)
; TVeg     = new((/case_sum,day_tot,ts_tot/),float,default_val)
; SMTemp   = new((/case_sum,day_tot,ts_tot/),float,default_val)

Tair     = 0.
Qle      = 0.
Qh       = 0.
deltaT   = 0.
Fwsoil   = 0.
Rnet     = 0.
Qg       = 0.
; TVeg     = 0.
; SMTemp   = 0.

do case_num = 0, case_sum -1
  print("case_num = " +case_num)

  do year = year_s,year_e
    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc" ; outputs/
    f1 = addfile (filename1,"r")

    Tair_tmp    = f1->Tair(0:60*8,{A:B},{C:D})
    Tair_tmp    = Tair_tmp - 273.15
    Qle_tmp     = f1->Qle(0:60*8,{A:B},{C:D})
    Qh_tmp      = f1->Qh(0:60*8,{A:B},{C:D})
    deltaT_tmp  = f1->VegT(0:60*8,{A:B},{C:D})
    deltaT_tmp  = deltaT_tmp - f1->Tair(0:60*8,{A:B},{C:D})
    Fwsoil_tmp  = f1->Fwsoil(0:60*8,{A:B},{C:D})

    Rnet_tmp    = f1->Rnet(0:60*8,{A:B},{C:D})
    ; Rnet_tmp    = Rnet_tmp + f1->Qh(0:60*8,{A:B},{C:D})
    Qg_tmp      = f1->Qg(0:60*8,{A:B},{C:D})
    ; TVeg_tmp    = f1->TVeg(0:60*8,{A:B},{C:D})
    ; SMTemp_tmp  = f1->SMTemp(0:60*8,{A:B},{C:D})

    printVarSummary(Qle_tmp)
    delete(f1)
    delete(filename1)

    counter = 0
    do event_num = 0,event_sum-1
      do day = HW_day_s(event_num),HW_day_e(event_num)

        day_s = (day-1)*8    ; the first day's 0h, since CABLE output 00:00 as 03:00,
                             ; 0 is 00:00
        day_e = (day-1)*8+8  ; the last day's 24h

        ; non-trees
        tmp = where(PFT_mask_multi .eq. 2, Tair_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Tair(case_num,counter,:)   = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 2, Qle_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Qle(case_num,counter,:)    = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 2, Qh_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        plot = gsn_csm_contour(wks1,dim_avg_n_Wrap(tmp,0),res1)
        Qh(case_num,counter,:)     = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 2, deltaT_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        deltaT(case_num,counter,:) = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 2, Fwsoil_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Fwsoil(case_num,counter,:) = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 2, Rnet_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Rnet(case_num,counter,:) = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 2, Qg_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Qg(case_num,counter,:)     = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        ; trees
        tmp = where(PFT_mask_multi .eq. 1, Tair_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Tair(case_num+case_sum,counter,:)   = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 1, Qle_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Qle(case_num+case_sum,counter,:)    = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 1, Qh_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        plot = gsn_csm_contour(wks1,dim_avg_n_Wrap(tmp,0),res1)
        Qh(case_num+case_sum,counter,:)     = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 1, deltaT_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        deltaT(case_num+case_sum,counter,:) = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 1, Fwsoil_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Fwsoil(case_num+case_sum,counter,:) = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 1, Rnet_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Rnet(case_num+case_sum,counter,:) = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 1, Qg_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Qg(case_num+case_sum,counter,:)     = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        counter = counter + 1
      end do
    end do

  end do ; year
  delete(Tair_tmp)
  delete(Qle_tmp)
  delete(Qh_tmp)
  delete(deltaT_tmp)
  delete(Fwsoil_tmp)

  delete(Rnet_tmp)
  delete(Qg_tmp)
  ; delete(TVeg_tmp)
  ; delete(SMTemp_tmp)

end do ; case_num
print("prepare to plot")
delete(PFT_mask_multi)

; 0, 1, 2, 3 , 4 , 5,  6,  7,  8
; 0, 3, 6, 9, 12, 15, 18, 21, 24

Tair_avg     = dim_avg_n_Wrap(Tair,1)
; Tair_min_max = new((/case_sum,2,ts_tot/),float,default_val)
; Tair_min_max(:,0,:) = dim_min_n_Wrap(Tair,1)
; Tair_min_max(:,1,:) = dim_max_n_Wrap(Tair,1)

Qle_avg     = dim_avg_n_Wrap(Qle,1)
; Qle_min_max = new((/case_sum,2,ts_tot/),float,default_val)
; Qle_min_max(:,0,:) = dim_min_n_Wrap(Qle,1)
; Qle_min_max(:,1,:) = dim_max_n_Wrap(Qle,1)

Qh_avg     = dim_avg_n_Wrap(Qh,1)
; Qh_min_max = new((/case_sum,2,ts_tot/),float,default_val)
; Qh_min_max(:,0,:) = dim_min_n_Wrap(Qh,1)
; Qh_min_max(:,1,:) = dim_max_n_Wrap(Qh,1)

deltaT_avg     = dim_avg_n_Wrap(deltaT,1)
; deltaT_min_max = new((/case_sum,2,ts_tot/),float,default_val)
; deltaT_min_max(:,0,:) = dim_min_n_Wrap(deltaT,1)
; deltaT_min_max(:,1,:) = dim_max_n_Wrap(deltaT,1)

Rnet_avg     = dim_avg_n_Wrap(Rnet,1)
; Rnet_min_max = new((/case_sum,2,ts_tot/),float,default_val)
; Rnet_min_max(:,0,:) = dim_min_n_Wrap(Rnet,1)
; Rnet_min_max(:,1,:) = dim_max_n_Wrap(Rnet,1)

EF_avg     = dim_avg_n_Wrap(Qle/(Qh+Qle),1)
; EF_min_max = new((/case_sum,2,ts_tot/),float,default_val)
; EF_min_max(:,0,:) = dim_min_n_Wrap(Qle/(Qh+Qle),1)
; EF_min_max(:,1,:) = dim_max_n_Wrap(Qle/(Qh+Qle),1)

Qg_avg     = dim_avg_n_Wrap(Qg,1)
; Qg_min_max = new((/case_sum,2,ts_tot/),float,default_val)
; Qg_min_max(:,0,:) = dim_min_n_Wrap(Qg,1)
; Qg_min_max(:,1,:) = dim_max_n_Wrap(Qg,1)

Fwsoil_avg     = dim_avg_n_Wrap(Fwsoil,(/1,2/))
print(Fwsoil_avg)

;_____________________________ PLOT ___________________________________
pic = "./plots/"+pic_name
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")

Thicknesses  = (/ 1.5, 1.5, 1.5, 1.5, \
                  1.5, 1.5, 1.5, 1.5  /)
DashPatterns = (/ 1, 1, 1, 1, \
                  0, 0, 0, 0  /)
Colors       = (/ "blue3","brown2","springgreen4","darkorange", \
                  "blue3","brown2","springgreen4","darkorange" /)
Opacities    = (/ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1/)

res                      = True
res@gsnDraw              = False
res@gsnFrame             = False
res@gsnMaximize          = False
res@vpHeightF            = 0.4
res@vpWidthF             = 0.70
res@vpXF                 = 0.08
res@vpYF                 = 0.7
res@tiXAxisFontHeightF   = 0.02
res@tiYAxisFontHeightF   = 0.02
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
res@gsnPaperOrientation  = "landscape";"portrait";"landscape" ;"auto"

x                        = ispan(1,ts_tot,1)
res@tmXBMode             = "Explicit"
res@tmXBValues           = ispan(1,ts_tot,1)
res@tmXBLabels           = tostring(ispan(0,24,3))
res@trXMinF              = 1.
res@trXMaxF              = ts_tot

res1 = res

res@xyLineThicknesses   = Thicknesses         ; make second line thicker
res@xyDashPatterns      = DashPatterns
res@xyLineColors        = Colors
res@tiXAxisString       = "hr"

; resources for "left" variable
resL                        = res
resL@tiYAxisString          = "Tleaf - Tair (~S~o~N~C)"
resL@trYMinF                = -4.
resL@trYMaxF                = 10.

resR                        = res
resR@tiYAxisString          = "Tair (~S~o~N~C)"
resR@trYMinF                = 20.
resR@trYMaxF                = 43.
resR@xyMarkLineMode         = "Markers"
resR@xyMarkers              = (/6, 16/)
resR@xyMarkerColors         = (/"red","red"/) ; 3 different colors
resR@xyExplicitLegendLabels = (/"non-Tree Tair","Tree Tair"/)

resR@pmLegendDisplayMode    = "Always"
resR@lgPerimOn              = False      ; No legend perimeter.
resR@pmLegendParallelPosF   = 0.25
resR@pmLegendOrthogonalPosF = -1.05
resR@lgBoxMinorExtentF      = 0.3
resR@pmLegendWidthF         = 0.20
resR@pmLegendHeightF        = 0.07
resR@lgItemOrder            = (/ 1,0 /)   ; Reorder the legends


plots = new(3,graphic)

plots(0) = gsn_csm_xy2(wks,x,deltaT_avg,Tair_avg(::case_sum,:),resL,resR)

res@pmLegendDisplayMode    = "Always"
res@lgPerimOn              = False
res@pmLegendParallelPosF   = 0.25
res@pmLegendOrthogonalPosF = -1.2
res@lgBoxMinorExtentF      = 0.3
res@pmLegendWidthF         = 0.20
res@pmLegendHeightF        = 0.25
res@lgItemOrder            = (/ 7,6,5,4,3,2,1,0 /)   ; Reorder the legends
res@xyExplicitLegendLabels = case_brf

res@tiYAxisString  = "LH (W m~S~-2~N~)"
res@trYMinF        = -50.
res@trYMaxF        = 400.
plots(1)           = gsn_csm_xy(wks,x,Qle_avg,res)


res@pmLegendDisplayMode = "Never"
res@tiYAxisString  = "SH (W m~S~-2~N~)"
res@trYMinF        = -50.
res@trYMaxF        = 400.
plots(2)           = gsn_csm_xy(wks,x,Qh_avg,res)

txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.028

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

do count = 0,2
  print("count = "+count)
  txid1   = gsn_create_text(wks, order_name(count), txres)
  annoid1 = gsn_add_annotation(plots(count), txid1, amres)
end do

pres                    = True
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
;pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01

gsn_panel(wks,(/plots/),(/3,1/),pres)

delete(plots)


plots = new(2,graphic)

res@tiYAxisString  = "Rnet (W m~S~-2~N~)"
res@trYMinF        =   -200.
res@trYMaxF        =   700.
plots(0)           = gsn_csm_xy(wks,x,Rnet_avg,res)

res@tiYAxisString  = "Qg (W m~S~-2~N~)"
res@trYMinF        = -600.
res@trYMaxF        = 600.
plots(1)           = gsn_csm_xy(wks,x,Qg_avg,res)


do count = 0,1
  print("count = "+count)
  txid1   = gsn_create_text(wks, order_name(count), txres)
  annoid1 = gsn_add_annotation(plots(count), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/2,1/),pres)

; ===================== shadow ========================
; Rnet_plots   = new(4,graphic)
; Qle_plots    = new(4,graphic)
; deltaT_plots = new(4,graphic)
; Qh_plots     = new(4,graphic)
;
; do case_num = 0, case_sum -1
;   print("case_num = "+case_num)
;   res1@xyLineThicknesses   = Thicknesses(case_num)
;   res1@gsnXYFillColors     = Colors(case_num)
;   res1@gsnXYFillOpacities  = Opacities(case_num)
;   res1@xyLineColor         = -1
;
;   ; Rnet
;   res1@trYMinF             =   -100.
;   res1@trYMaxF             =   600.
;   Rnet_plots(case_num)   = gsn_csm_xy(wks,x,Rnet_min_max(case_num,:,:),res1)
;
;   ; Qle
;   res1@trYMinF             = -100.
;   res1@trYMaxF             = 600.
;   Qle_plots(case_num)  = gsn_csm_xy(wks,x,Qle_min_max(case_num,:,:),res1)
;
;   ; deltaT
;   res1@trYMinF             = -5.
;   res1@trYMaxF             = 10.
;   deltaT_plots(case_num)  = gsn_csm_xy(wks,x,deltaT_min_max(case_num,:,:),res1)
;
;   ; Qh
;   res1@trYMinF              = -100.
;   res1@trYMaxF              = 600.
;   Qh_plots(case_num)  = gsn_csm_xy(wks,x,Qh_min_max(case_num,:,:),res1)
;
;   if case_num .eq. 0 then
;     overlay(Rnet_plots(case_num),   plots(0))
;     overlay(Qle_plots(case_num),    plots(1))
;     overlay(deltaT_plots(case_num), plots(2))
;     overlay(Qh_plots(case_num),     plots(3))
;   else
;     overlay(Rnet_plots(case_num),   Rnet_plots(case_num-1))
;     overlay(Qle_plots(case_num),    Qle_plots(case_num-1))
;     overlay(deltaT_plots(case_num), deltaT_plots(case_num-1))
;     overlay(Qh_plots(case_num),     Qh_plots(case_num-1))
;   end if
; end do
; plots1(0) = Rnet_plots(3)
; plots1(1) = Qle_plots(3)
; plots1(2) = deltaT_plots(3)
; plots1(3) = Qh_plots(3)


end
