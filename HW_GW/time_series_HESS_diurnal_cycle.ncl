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
pic_name   = "CABLE-HESS_"+message+"_PFT-"+PFT+"_"+scale+"_"+year_s+"-"+year_e

path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"

case_name = (/"CTL-30x3+30yr","FREE_DRAIN", "CTL-30x3+30yr_2019-WetSoil",\
              "CTL-30x3+30yr_2019-DeepRoot"/)

order_name  = (/"(a)","(b)","(c)","(d)","(e)","(f)",\
                "(g)","(h)","(i)","(j)","(k)","(l)"/)

; 0, 1, 2, 3 , 4 , 5,  6,  7
; 3, 6, 9, 12, 15, 18, 21, 24

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

;______________________________ Data _________________________________
case_sum  = dimsizes(case_name)
path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = path_mum+case_name(case_num)+"/outputs-subday"
end do

; ___________________________ TREE Mask ______________________________
filename1 = path(0)+"/cable_out_2019_SE_Aus.nc"
f1 = addfile (filename1,"r")
PFT_mask = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
default_val = f1->Qle@_FillValue

PFT_mask_multi = new((/ts_tot,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
do z = 0, ts_tot-1
  PFT_mask_multi(z,:,:) = PFT_mask
end do

delete(filename1)
delete(f1)

Tair     = new((/case_sum,day_tot,ts_tot/),float,default_val)
Qle      = new((/case_sum,day_tot,ts_tot/),float,default_val)
Qh       = new((/case_sum,day_tot,ts_tot/),float,default_val)
deltaT   = new((/case_sum,day_tot,ts_tot/),float,default_val)
Fwsoil   = new((/case_sum,day_tot,ts_tot/),float,default_val)
Rnet     = new((/case_sum,day_tot,ts_tot/),float,default_val)
; Qg       = new((/case_sum,day_tot,ts_tot/),float,default_val)
; TVeg     = new((/case_sum,day_tot,ts_tot/),float,default_val)
; SMTemp   = new((/case_sum,day_tot,ts_tot/),float,default_val)

Tair     = 0.
Qle      = 0.
Qh       = 0.
deltaT   = 0.
Fwsoil   = 0.
Rnet     = 0.
; Qg       = 0.
; TVeg     = 0.
; SMTemp   = 0.

do case_num = 0, case_sum -1
  print("case_num = " +case_num)

  do year = year_s,year_e
    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc" ; outputs/
    f1 = addfile (filename1,"r")

    Tair_tmp    = f1->Tair(0:60*8,{A:B},{C:D})
    Qle_tmp     = f1->Qle(0:60*8,{A:B},{C:D})
    Qh_tmp      = f1->Qh(0:60*8,{A:B},{C:D})
    deltaT_tmp  = f1->VegT(0:60*8,{A:B},{C:D})
    deltaT_tmp  = deltaT_tmp - Tair_tmp
    Fwsoil_tmp  = f1->Fwsoil(0:60*8,{A:B},{C:D})

    Rnet_tmp    = f1->Qle(0:60*8,{A:B},{C:D})
    Rnet_tmp    = Rnet_tmp + f1->Qh(0:60*8,{A:B},{C:D})
    ; Qg_tmp      = f1->Qg(0:60*8,{A:B},{C:D})
    ; TVeg_tmp    = f1->TVeg(0:60*8,{A:B},{C:D})
    ; SMTemp_tmp  = f1->SMTemp(0:60*8,{A:B},{C:D})

    printVarSummary(Qle_tmp)
    delete(f1)
    delete(filename1)

    counter = 0
    do event_num = 0,event_sum-1
      do day = HW_day_s(event_num),HW_day_e(event_num)

        day_s = (day-1)*8-1 ; the first day's 0h
        day_e = (day-1)*8+7 ; the last day's 24h

        tmp = where(PFT_mask_multi .eq. 1, Tair_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Tair(case_num,counter,:)   = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 1, Qle_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Qle(case_num,counter,:)    = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 1, Qh_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Qh(case_num,counter,:)     = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 1, deltaT_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        deltaT(case_num,counter,:) = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 1, Fwsoil_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Fwsoil(case_num,counter,:) = dim_avg_n_Wrap(tmp,(/1,2/))
        delete(tmp)

        tmp = where(PFT_mask_multi .eq. 1, Rnet_tmp(day_s:day_e,:,:), default_val)
        tmp@_FillValue = default_val
        Rnet(case_num,counter,:) = dim_avg_n_Wrap(tmp,(/1,2/))
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
  ; delete(Qg_tmp)
  ; delete(TVeg_tmp)
  ; delete(SMTemp_tmp)

end do ; case_num
print("prepare to plot")
delete(PFT_mask_multi)

; 0, 1, 2, 3 , 4 , 5,  6,  7,  8
; 0, 3, 6, 9, 12, 15, 18, 21, 24

Tair_avg     = dim_avg_n_Wrap(Tair,1)
Tair_min_max = new((/case_sum,2,ts_tot/),float,default_val)
Tair_min_max(:,0,:) = dim_min_n_Wrap(Tair,1)
Tair_min_max(:,1,:) = dim_max_n_Wrap(Tair,1)

Qle_avg     = dim_avg_n_Wrap(Qle,1)
Qle_min_max = new((/case_sum,2,ts_tot/),float,default_val)
Qle_min_max(:,0,:) = dim_min_n_Wrap(Qle,1)
Qle_min_max(:,1,:) = dim_max_n_Wrap(Qle,1)

Qh_avg     = dim_avg_n_Wrap(Qh,1)
Qh_min_max = new((/case_sum,2,ts_tot/),float,default_val)
Qh_min_max(:,0,:) = dim_min_n_Wrap(Qh,1)
Qh_min_max(:,1,:) = dim_max_n_Wrap(Qh,1)

deltaT_avg     = dim_avg_n_Wrap(deltaT,1)
deltaT_min_max = new((/case_sum,2,ts_tot/),float,default_val)
deltaT_min_max(:,0,:) = dim_min_n_Wrap(deltaT,1)
deltaT_min_max(:,1,:) = dim_max_n_Wrap(deltaT,1)

Rnet_avg     = dim_avg_n_Wrap(Rnet,1)
Rnet_min_max = new((/case_sum,2,ts_tot/),float,default_val)
Rnet_min_max(:,0,:) = dim_min_n_Wrap(Rnet,1)
Rnet_min_max(:,1,:) = dim_max_n_Wrap(Rnet,1)

Fwsoil_avg     = dim_avg_n_Wrap(Fwsoil,(/1,2/))


;_____________________________ PLOT ___________________________________
pic = "./plots/"+pic_name
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

res                    = True
res@gsnDraw            = True
res@gsnFrame           = True
res@gsnMaximize        = True
res@vpHeightF          = 0.5
res@vpWidthF           = 0.90
res@vpXF               = 0.08
res@vpYF               = 0.7
res@pmLegendWidthF     = 0.10                     ; Resize legend width
res@pmLegendHeightF    = 0.10                     ; and height
res@tiXAxisFontHeightF = 0.02
res@tiYAxisFontHeightF = 0.02
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
res@gsnPaperOrientation = "landscape";"portrait";"landscape" ;"auto"

Thicknesses = (/  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0/)
DashPatterns= (/0, 0, 0, 0, 0, 0, 0, 0, 0/)
Colors      = (/"black","red","blue","green","purple","blue","chartreuse","orange","cyan","yellow"/)
Opacities   = (/0.1, 0.1, 0.1, 0.1/)



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

res@tmXBMode  = "Explicit"
x             = ispan(1,ts_tot,1)
res@tmXBValues= ispan(1,ts_tot,1)
res@tmXBLabels= tostring(ispan(0,24,3))
res@trXMinF   = 1.
res@trXMaxF   = ts_tot

res1 = res

res@xyLineThicknesses   = Thicknesses         ; make second line thicker
res@xyDashPatterns      = DashPatterns
res@xyLineColors        = Colors
res@xyExplicitLegendLabels = case_name

res@tiXAxisString  = "hr"
plots = new(4,graphic)

res@pmLegendDisplayMode = "Never"

res@tiYAxisString  = "Rnet (W/m2)"
res@trYMinF        =   -100.
res@trYMaxF        =   600.
plots(0)           = gsn_csm_xy(wks,x,Rnet_avg,res)

res@tiYAxisString  = "Qle (W/m2)"
res@trYMinF        = -100.
res@trYMaxF        = 600.
plots(1)           = gsn_csm_xy(wks,x,Qle_avg,res)

res@tiYAxisString  = "DeltaT (C)"
res@trYMinF        = -20.
res@trYMaxF        = 20.
plots(2)           = gsn_csm_xy(wks,x,deltaT_avg,res)

res@tiYAxisString  = "Qh (W/m2)"
res@trYMinF        = -100.
res@trYMaxF        = 600.
plots(3)           = gsn_csm_xy(wks,x,Qh_avg,res)

Rnet_plots   = new(4,graphic)
Qle_plots    = new(4,graphic)
deltaT_plots = new(4,graphic)
Qh_plots     = new(4,graphic)

do case_num = 0, case_sum -1
  print("case_num = "+case_num)
  res1@xyLineThicknesses   = Thicknesses(case_num)
  res1@gsnXYFillColors     = Colors(case_num)
  res1@gsnXYFillOpacities  = Opacities(case_num)
  res1@xyLineColor         = -1

  ; Rnet
  res1@trYMinF             =   -100.
  res1@trYMaxF             =   600.
  Rnet_plots(case_num)   = gsn_csm_xy(wks,x,Rnet_min_max(case_num,:,:),res1)

  ; Qle
  res1@trYMinF             = -100.
  res1@trYMaxF             = 600.
  Qle_plots(case_num)  = gsn_csm_xy(wks,x,Qle_min_max(case_num,:,:),res1)

  ; deltaT
  res1@trYMinF             = -20.
  res1@trYMaxF             = 20.
  deltaT_plots(case_num)  = gsn_csm_xy(wks,x,deltaT_min_max(case_num,:,:),res1)

  ; Qh
  res1@trYMinF              = -100.
  res1@trYMaxF              = 600.
  Qh_plots(case_num)  = gsn_csm_xy(wks,x,Qh_min_max(case_num,:,:),res1)

  if case_num .eq. 0 then
    overlay(Rnet_plots(case_num),   plots(0))
    overlay(Qle_plots(case_num),    plots(1))
    overlay(deltaT_plots(case_num), plots(2))
    overlay(Qh_plots(case_num),     plots(3))
  else
    overlay(Rnet_plots(case_num),   Rnet_plots(case_num-1))
    overlay(Qle_plots(case_num),    Qle_plots(case_num-1))
    overlay(deltaT_plots(case_num), deltaT_plots(case_num-1))
    overlay(Qh_plots(case_num),     Qh_plots(case_num-1))    
  end if
end do

txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.03

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

do count = 0,3
  txid1   = gsn_create_text(wks, order_name(count), txres)
  annoid1 = gsn_add_annotation(plots(count), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/4,1/),pres)
end

;
; res1 = res
; res1@tiYAxisString  = "Rain (mm/3hr)"
; res1@gsnXYBarChart  = True
; ; res1@tmYUseLeft     = False
; ; res1@tmYLOn         = False
; ; res1@tmYROn         = True
; ; res1@tmYRLabelsOn   = True
; res1@trYMinF        = 0.
; res1@trYMaxF        = 2.
; res1@xyLineColor    = "gray"
; ; plots(0) = gsn_csm_xy(wks,x,Rain(0,:),res1)

; res@tiYAxisString  = "SMTemp-Tair (C)"
; res@trYMinF        =   0.
; res@trYMaxF        =  10.2
; plots(1)           = gsn_csm_xy(wks,x,SMTemp-Tair,res)
;
; res@tiYAxisString  = "Rnet (W/m2)"
; res@trYMinF        = -100.
; res@trYMaxF        = 700.
; plots(2)           = gsn_csm_xy(wks,x,Rnet,res)
;
; res@tiYAxisString  = "Qg (W/m2)"
; res@trYMinF        = -400.
; res@trYMaxF        = 400.
; plots(5)           = gsn_csm_xy(wks,x,Qg,res)

; res@tiYAxisString  = "deltaT (C)"
; res@trYMinF        =   0.
; res@trYMaxF        =  10.2
; plots(6) = gsn_csm_xy(wks,x,deltaT,res)
;
; res@pmLegendDisplayMode = "Always"
; res@tiYAxisString  = "TVeg (mm/3hr)"
; res@trYMinF        =   0.
; res@trYMaxF        =   1.
; plots(7) = gsn_csm_xy(wks,x,TVeg,res)
;