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
time_scale = "daily" ;"monthly" ; "daily"
scale      = "SE-AU"
year_s     = 2000
year_e     = 2019
year_sum   = year_e-year_s+1

; ______________________________ Path _________________________________
path_mum       = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
case_name      = (/"CTL-30x3+30yr/outputs", "FREE_DRAIN/outputs-spinup30x3yr"/) ;,"CTL-30x3+30yr-DeepRoot/outputs"/)
pic_message    = "CTL_vs_FREE_DRAIN"
pic_name       = "HESS_EF_"+pic_message+"_"+time_scale+"_"+scale+"_"+year_s+"-"+year_e
Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Depth          = (/1.1, 5.1, 15.7, 43.85, 118.55, 316.4, 2355.36/)
                 ; avg(GWdz) over SE Aus = 37.9072 m
                 ; -> last depth=460+(37.9072*10)/2 = 2355.36 (cm)

order_name     = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)

; "tree" PFT
iveg_s         = 1
iveg_e         = 4

if time_scale .eq. "monthly" then
  print("monthly")
  time_tot = year_sum*12
else if time_scale .eq. "daily" then
  print("daily")
  time_tot = 0
  do year = year_s,year_e
    if mod(year,4) .eq. 0 then
      time_tot = time_tot + 366
    else
      time_tot = time_tot + 365
    end if
  end do
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
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
    f1 = addfile (filename1,"r")

    if time_scale .eq. "monthly" then
      ts = 12
    else if time_scale .eq. "daily" then
      ts = sum(dom)
    end if
    end if

    PFT_mask_multi = new((/ts,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
    do z = 0, ts-1
      PFT_mask_multi(z,:,:) = PFT_mask
    end do

    Qle(case_num, 0, i:i+ts-1)  = dim_avg_n_Wrap(f1->Qle(:,{A:B},{C:D})*1., (/1,2/))
    Qh(case_num, 0, i:i+ts-1)   = dim_avg_n_Wrap(f1->Qh(:,{A:B},{C:D})*1., (/1,2/))
    Rnet(case_num, 0, i:i+ts-1) = dim_avg_n_Wrap(f1->Rnet(:,{A:B},{C:D})*1., (/1,2/))
    Qle(case_num, 1, i:i+ts-1)  = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qle(:,{A:B},{C:D}),default_val), (/1,2/))
    Qh(case_num, 1, i:i+ts-1)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qh(:,{A:B},{C:D}),default_val), (/1,2/))
    Rnet(case_num, 1, i:i+ts-1) = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Rnet(:,{A:B},{C:D}),default_val), (/1,2/))

    i = i + ts
    delete(PFT_mask_multi)
  end do ; year
end do ; case_num
print("prepare to plot")

EF_all  = new((/case_sum*2,time_tot/),float, default_val)
EF_tree = new((/case_sum*2,time_tot/),float, default_val)
Rnet_thd = 0.
EF_all(0,:)  = where( Rnet(0,0,:) .gt. Rnet_thd, Qle(0,0,:)/Rnet(0,0,:), default_val) ; CTL
EF_all(1,:)  = where( Rnet(1,0,:) .gt. Rnet_thd, Qle(1,0,:)/Rnet(1,0,:), default_val) ; FREE_DRAIN
EF_all(2,:)  = where( Rnet(0,0,:) .gt. Rnet_thd, Qh(0,0,:)/Rnet(0,0,:), default_val)
EF_all(3,:)  = where( Rnet(1,0,:) .gt. Rnet_thd, Qh(1,0,:)/Rnet(1,0,:), default_val)

EF_tree(0,:)  = where( Rnet(0,1,:) .gt. Rnet_thd, Qle(0,1,:)/Rnet(0,1,:), default_val) ; CTL
EF_tree(1,:)  = where( Rnet(1,1,:) .gt. Rnet_thd, Qle(1,1,:)/Rnet(1,1,:), default_val) ; FREE_DRAIN
EF_tree(2,:)  = where( Rnet(0,1,:) .gt. Rnet_thd, Qh(0,1,:)/Rnet(0,1,:), default_val)
EF_tree(3,:)  = where( Rnet(1,1,:) .gt. Rnet_thd, Qh(1,1,:)/Rnet(1,1,:), default_val)
EF_all  = runave(EF_all,90,0)
EF_tree = runave(EF_tree,90,0)
delete(Qle)
delete(Qh)

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

res@xyLineThicknesses       = (/ 1.5, 1.5, 1.5, 1.5/)
res@xyDashPatterns          = (/   0,   0,   0,   0/)
res@xyLineColors            = (/"blue","green","red","orange"/)
res@xyExplicitLegendLabels  = (/"CTL_Qle/Rnet","FD_Qle/Rnet","CTL_Qh/Rnet","FD_Qh/Rnet"/)

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

if time_scale .eq. "monthly" then
  res@tmXBValues            = ispan(1,time_tot,12*2)
  res@tmXBLabels            = ispan(year_s,year_e,2)
  res@tiXAxisString         = "Year"
else if time_scale .eq. "daily" then
  res@tmXBValues            = ispan(1,time_tot,366*2)
  res@tmXBLabels            = ispan(year_s,year_e,2)
  res@tiXAxisString         = "Year"
end if
end if

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


; _______________________________ plot soil profile _____________________________
plots              = new(2,graphic)
x                  = ispan(1,time_tot,1)
res@trYMinF        = 0.
res@trYMaxF        = 1.1
res@pmLegendDisplayMode = "Always"
plots(0)           = gsn_csm_xy(wks, x, EF_all,res)
res@pmLegendDisplayMode = "Never"
plots(1)           = gsn_csm_xy(wks, x, EF_tree,res)

do counter = 0,1
  txid1   = gsn_create_text(wks, order_name(counter), txres)
  annoid1 = gsn_add_annotation(plots(counter), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/1,2/),pres)

end
