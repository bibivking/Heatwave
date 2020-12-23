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
HW_day_s   = 11
HW_day_e   = 11
year_sum   = year_e-year_s+1
time_tot   = year_sum*365*8
plot_ts_s  = (HW_day_s-1)*8-1 ; the first day's 0h
plot_ts_e  = (HW_day_e-1)*8+7 ; the last day's 24h
plot_ts    = plot_ts_e - plot_ts_s + 1

message    = "subday-24h"
pic_name   = "CABLE-HESS_"+message+"_PFT-"+PFT+"_"+scale+"_"+year_s+"-"+year_e

path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"

case_name = (/"CTL-30x3+30yr","FREE_DRAIN", "CTL-30x3+30yr_2019-WetSoil",\
              "CTL-30x3+30yr_2019-DeepRoot"/)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers


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

PFT_mask_multi = new((/plot_ts,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
do z = 0, plot_ts-1
  PFT_mask_multi(z,:,:) = PFT_mask
end do

delete(filename1)
delete(f1)

Rnet     = new((/case_sum,plot_ts/),float,default_val)
Rain     = new((/case_sum,plot_ts/),float,default_val)
Tair     = new((/case_sum,plot_ts/),float,default_val)
Qle      = new((/case_sum,plot_ts/),float,default_val)
Qh       = new((/case_sum,plot_ts/),float,default_val)
Qg       = new((/case_sum,plot_ts/),float,default_val)
TVeg     = new((/case_sum,plot_ts/),float,default_val)
SMTemp   = new((/case_sum,plot_ts/),float,default_val)
deltaT   = new((/case_sum,plot_ts/),float,default_val)
Fwsoil   = new((/case_sum,plot_ts/),float,default_val)

Rnet     = 0.
Rain     = 0.
Tair     = 0.
Qle      = 0.
Qh       = 0.
Qg       = 0.
SMTemp   = 0.
TVeg     = 0.
deltaT   = 0.
Fwsoil   = 0.

do case_num = 0, case_sum -1
  print("case_num = " +case_num)
  do year = year_s,year_e
    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc" ; outputs/
    f1 = addfile (filename1,"r")
    print(f1->time(plot_ts_s)/3600)
    print(f1->time(plot_ts_e)/3600)
    ; tmp = f1->Qle(plot_ts_s:plot_ts_e,{A:B},{C:D})/(f1->Qle(plot_ts_s:plot_ts_e,{A:B},{C:D})+f1->Qh(plot_ts_s:plot_ts_e,{A:B},{C:D}))
    ; Qle(case_num,:)      = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1 .and. tmp .gt. 0, tmp, default_val) ,(/1,2/))
    Rnet(case_num,:)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Rnet(plot_ts_s:plot_ts_e,{A:B},{C:D}),  default_val) ,(/1,2/))
    Rain(case_num,:)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Rainf(plot_ts_s:plot_ts_e,{A:B},{C:D}), default_val) ,(/1,2/))*3600*3
    Tair(case_num,:)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Tair(plot_ts_s:plot_ts_e,{A:B},{C:D}),  default_val) ,(/1,2/))-273.15
    Qle(case_num,:)      = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qle(plot_ts_s:plot_ts_e,{A:B},{C:D}),   default_val) ,(/1,2/))
    Qh(case_num,:)       = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qh(plot_ts_s:plot_ts_e,{A:B},{C:D}),    default_val) ,(/1,2/))
    Qg(case_num,:)       = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qg(plot_ts_s:plot_ts_e,{A:B},{C:D}),    default_val) ,(/1,2/))
    TVeg(case_num,:)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->TVeg(plot_ts_s:plot_ts_e,{A:B},{C:D}),  default_val) ,(/1,2/))*3600*3
    SMTemp(case_num,:)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->SoilTemp(plot_ts_s:plot_ts_e,0,{A:B},{C:D}),  default_val) ,(/1,2/)) - 273.15
    deltaT(case_num,:)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->VegT(plot_ts_s:plot_ts_e,{A:B},{C:D}) - f1->Tair(plot_ts_s:plot_ts_e,{A:B},{C:D}), default_val) ,(/1,2/))
    ; deltaT(case_num*3,:)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->VegT(plot_ts_s:plot_ts_e,{A:B},{C:D}),  default_val) ,(/1,2/)) - 273.15
    ; deltaT(case_num*3+1,:) = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Tair(plot_ts_s:plot_ts_e,{A:B},{C:D}),  default_val) ,(/1,2/)) - 273.15
    ; deltaT(case_num*3+2,:) = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->VegT(plot_ts_s:plot_ts_e,{A:B},{C:D})-f1->Tair(plot_ts_s:plot_ts_e,{A:B},{C:D}),  default_val) ,(/1,2/))
    Fwsoil(case_num,:)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Fwsoil(plot_ts_s:plot_ts_e,{A:B},{C:D}),default_val) ,(/1,2/))
  end do ; year
  delete(filename1)
  delete(f1)
end do ; case_num
print("prepare to plot")
delete(PFT_mask_multi)



; 0, 1, 2, 3 , 4 , 5,  6,  7,  8
; 0, 3, 6, 9, 12, 15, 18, 21, 24
count = 0
do day = HW_day_s,HW_day_e
  ; Rain(:,count*8:count*8+1) = default_val
  ; Rnet(:,count*8+1:count*8+2)   = default_val
  ; Tair(:,count*8+1:count*8+2)   = default_val
  ; Qle(:,count*8+1:count*8+2)    = default_val
  ; Qh(:,count*8+1:count*8+2)     = default_val
  ; Qg(:,count*8+1:count*8+2)     = default_val
  ; TVeg(:,count*8+1:count*8+2)   = default_val
  ; deltaT(:,count*8+1:count*8+2) = default_val
  ; Fwsoil(:,count*8+1:count*8+2) = default_val

  ; Rain(:,count*8+6:count*8+7)  = default_val
  ; Tair(:,count*8+6:count*8+7)  = default_val
  ; Qle(:,count*8+7)   = default_val
  ; Qh(:,count*8+7)    = default_val
  ; TVeg(:,count*8+7)  = default_val
  ; deltaT(:,count*8+7)= default_val
  ; Fwsoil(:,count*8+7)= default_val
  count = count + 1
end do

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

res@xyLineThicknesses   = (/  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0,1.0/)          ; make second line thicker
res@xyDashPatterns      = (/0, 0, 0, 0, 0, 0, 0, 0, 0/)
res@xyLineColors        = (/"black","red","blue","green","purple","blue","chartreuse","orange","cyan","yellow"/)          ; change line color
res@xyExplicitLegendLabels = case_name
;res@trYMinF = 0

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


; Qle_diff    = new((/case_sum-1,plot_ts/),float,default_val)
; Qh_diff     = new((/case_sum-1,plot_ts/),float,default_val)
; TVeg_diff   = new((/case_sum-1,plot_ts/),float,default_val)
; deltaT_diff = new((/case_sum-1,plot_ts/),float,default_val)
; Fwsoil_diff = new((/case_sum-1,plot_ts/),float,default_val)
;
; Qle_diff(0,:) =  Qle(1,:) - Qle(0,:)
; Qle_diff(1,:) =  Qle(2,:) - Qle(0,:)
; Qle_diff(2,:) =  Qle(3,:) - Qle(0,:)
;
; Qh_diff(0,:) =  Qh(1,:) - Qh(0,:)
; Qh_diff(1,:) =  Qh(2,:) - Qh(0,:)
; Qh_diff(2,:) =  Qh(3,:) - Qh(0,:)
;
; TVeg_diff(0,:) =  TVeg(1,:) - TVeg(0,:)
; TVeg_diff(1,:) =  TVeg(2,:) - TVeg(0,:)
; TVeg_diff(2,:) =  TVeg(3,:) - TVeg(0,:)
;
; deltaT_diff(0,:) =  deltaT(1,:) - deltaT(0,:)
; deltaT_diff(1,:) =  deltaT(2,:) - deltaT(0,:)
; deltaT_diff(2,:) =  deltaT(3,:) - deltaT(0,:)
;
; Fwsoil_diff(0,:) =  Fwsoil(1,:) - Fwsoil(0,:)
; Fwsoil_diff(1,:) =  Fwsoil(2,:) - Fwsoil(0,:)
; Fwsoil_diff(2,:) =  Fwsoil(3,:) - Fwsoil(0,:)

res@tmXBMode  = "Explicit"
x             = ispan(1,plot_ts,1)
res@tmXBValues= ispan(1,plot_ts,1)
res@tmXBLabels= tostring(ispan(0,24,3))

; res@tmXBValues= ispan(1,plot_ts,8)
; res@tmXBLabels= tostring(ispan(HW_day_s,HW_day_e+1,1))
res@tiXAxisString  = "day"
plots = new(8,graphic)

res@pmLegendDisplayMode = "Never"

res1 = res
res1@tiYAxisString  = "Rain (mm/3hr)"
res1@gsnXYBarChart  = True
; res1@tmYUseLeft     = False
; res1@tmYLOn         = False
; res1@tmYROn         = True
; res1@tmYRLabelsOn   = True
res1@trYMinF        = 0.
res1@trYMaxF        = 2.
res1@xyLineColor    = "gray"
; plots(0) = gsn_csm_xy(wks,x,Rain(0,:),res1)

res@tiYAxisString  = "Tair (C)"
res@trYMinF        = 15.
res@trYMaxF        = 45.
plots(0)           = gsn_csm_xy(wks,x,Tair(0,:),res)

res@tiYAxisString  = "SMTemp-Tair (C)"
res@trYMinF        =   0.
res@trYMaxF        =  10.2
plots(1)           = gsn_csm_xy(wks,x,SMTemp-Tair,res)

res@tiYAxisString  = "Rnet (W/m2)"
res@trYMinF        = -100.
res@trYMaxF        = 700.
plots(2)           = gsn_csm_xy(wks,x,Rnet,res)
; overlay(plots(0),plot)

res@tiYAxisString  = "Qle (W/m2)"
res@trYMinF        =   -100.
res@trYMaxF        = 400.
plots(3) = gsn_csm_xy(wks,x,Qle,res)

res@tiYAxisString  = "Qh (W/m2)"
res@trYMinF        = -100.
res@trYMaxF        = 400.
plots(4) = gsn_csm_xy(wks,x,Qh,res)


res@tiYAxisString  = "Qg (W/m2)"
res@trYMinF        = -400.
res@trYMaxF        = 400.
plots(5)           = gsn_csm_xy(wks,x,Qg,res)

res@tiYAxisString  = "deltaT (C)"
res@trYMinF        =   0.
res@trYMaxF        =  10.2
plots(6) = gsn_csm_xy(wks,x,deltaT,res)

res@pmLegendDisplayMode = "Always"
res@tiYAxisString  = "TVeg (mm/3hr)"
res@trYMinF        =   0.
res@trYMaxF        =   1.
plots(7) = gsn_csm_xy(wks,x,TVeg,res)

; res@tiYAxisString  = "Fwsoil (-)"
; res@trYMinF        =   0.
; res@trYMaxF        =   1.
; plots(5) = gsn_csm_xy(wks,x,Fwsoil,res)

txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.03

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

do count = 0,5
  txid1   = gsn_create_text(wks, order_name(count), txres)
  annoid1 = gsn_add_annotation(plots(count), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/4,2/),pres)
end
