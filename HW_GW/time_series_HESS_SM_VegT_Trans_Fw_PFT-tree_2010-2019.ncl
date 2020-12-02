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
PFT        = "tree" ; "tree" ;"all"
year_s     = 2010
year_e     = 2019
year_sum   = year_e-year_s+1

; ______________________________ Path _________________________________
path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
case_name = (/"CTL","FREE_DRAIN"/);
           ;"FREE_DRAIN","PUMPx10","PUMPx100"/)
           ;,"PUMP","PUMP_IRRG","PUMP_IRRGx10","PUMP_IRRGx100"/)
pic_message = "no-spinup"
pic_name   = "HESS_SM_VegT_Trans_Fw_"+case_name(0)+"_vs_"+case_name(1)+"_"+pic_message+"_PFT-"+PFT+"_"+time_scale+"_"+scale+"_"+year_s+"-"+year_e
Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Depth          = (/1.1, 5.1, 15.7, 43.85, 118.55, 316.4, 1710/) ; including aquifer

if time_scale .eq. "monthly" then
  print("monthly")
  time_tot = year_sum*12
else if time_scale .eq. "daily" then
  print("daily")
  time_tot = year_sum*365 + 2
  print(time_tot)
end if
end if

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
  path(case_num) = path_mum+case_name(case_num)+"/outputs-without-spinup"
end do

; ___________________________ TREE Mask ______________________________
filename1 = path(0)+"/cable_out_1990_SE_Aus.nc"
f1 = addfile (filename1,"r")
PFT_mask = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)

default_val = f1->Qle@_FillValue

delete(filename1)
delete(f1)

; Rain     = new((/case_sum,time_tot/),float,default_val)
SoilMoist = new((/case_sum,7,time_tot/),float,default_val)
SM64cm    = new((/case_sum,time_tot/),float,default_val)
Qrec      = new((/case_sum+1,time_tot/),float,default_val)
VegT     = new((/case_sum,time_tot/),float,default_val)
TVeg     = new((/case_sum,time_tot/),float,default_val)
Fwsoil   = new((/case_sum,time_tot/),float,default_val)
; EVAP64cm = new((/case_sum,time_tot/),float,default_val)
; EVAPFBL  = new((/case_sum,6,time_tot/),float,default_val)

; Rain     = 0.
SoilMoist= 0.
SM64cm   = 0.
Qrec     = 0.
VegT     = 0.
TVeg     = 0.
Fwsoil   = 0.
; EVAP64cm = 0.
; EVAPFBL  = 0.

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

    ; Rain(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Rainf(:,{A:B},{C:D}),default_val) ,(/1,2/))*3600*24
    VegT(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->VegT(:,{A:B},{C:D}),default_val) ,(/1,2/)) -273.15
    TVeg(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->TVeg(:,{A:B},{C:D}),default_val) ,(/1,2/))*3600*24
    Fwsoil(case_num,i:i+ts-1)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Fwsoil(:,{A:B},{C:D}),default_val) ,(/1,2/))

    tmp = ( f1->SoilMoist(:,0,{A:B},{C:D})*0.022 + f1->SoilMoist(:,1,{A:B},{C:D})*0.058 \
          + f1->SoilMoist(:,2,{A:B},{C:D})*0.154 + f1->SoilMoist(:,3,{A:B},{C:D})*0.409)/0.643
    SM64cm(case_num,i:i+ts-1) = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, tmp, default_val), (/1,2/))

    ; tmp1 = f1->EVAPFBL(:,0,{A:B},{C:D}) + f1->EVAPFBL(:,1,{A:B},{C:D}) \
    ;       + f1->EVAPFBL(:,2,{A:B},{C:D}) + f1->EVAPFBL(:,3,{A:B},{C:D})
    ; EVAP64cm(case_num,i:i+ts-1) = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, tmp1*3600*24, default_val), (/1,2/))

    Qrec(case_num,i:i+ts-1)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qrecharge(:,{A:B},{C:D}),default_val) ,(/1,2/))*3600*24
    do soil = 0,5
       SoilMoist(case_num, soil, i:i+ts-1)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->SoilMoist(:,soil,{A:B},{C:D}),default_val), (/1,2/))
       ; EVAPFBL(case_num, soil, i:i+ts-1)  = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->EVAPFBL(:,soil,{A:B},{C:D})*3600*24,default_val), (/1,2/))
    end do
    ; aquifer
    SoilMoist(case_num, 6, i:i+ts-1)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->GWMoist(:,{A:B},{C:D}),default_val), (/1,2/))
    i = i + ts
    delete(tmp)
    ; delete(tmp1)
    delete(PFT_mask_multi)
  end do ; year
end do ; case_num
print("prepare to plot")

Qrec = where(Qrec .lt. 0, 1, -1)

;_____________________________ PLOT ___________________________________
pic = "./plots/"+pic_name
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��
res                    = True
res@gsnMaximize        = True
res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ����ҳ
res@vpHeightF          = 0.4
res@vpWidthF           = 0.88
res@vpXF               = 0.08
res@vpYF               = 0.7
res@pmLegendWidthF     = 0.2                     ; Resize legend width
res@pmLegendHeightF    = 0.2                     ; and height
res@tiXAxisFontHeightF = 0.02
res@tiYAxisFontHeightF = 0.02
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
res@gsnPaperOrientation = "landscape";"portrait";"landscape" ;"auto"

res@tmXBMode           = "Explicit"
res@tmBorderThicknessF = 3.0
if time_scale .eq. "monthly" then
  res@tmXBValues= ispan(1,time_tot,12)
  res@tiXAxisString  = "months"
else if time_scale .eq. "daily" then
  res@tmXBValues= ispan(1,time_tot,366)
  res@tiXAxisString  = "days"
end if
end if

res@tmXBLabels= tostring(ispan(year_s,year_e,1))

res1                    = res
res3                    = res

res@xyLineThicknesses   = (/  1.5, 1.5, 1.5, 1.5, 1.5, 1.5/) ;, 1.5, 1.5, 1.5,1.5/)          ; make second line thicker
res@pmLegendDisplayMode = "never" ;"Always"
res@xyDashPatterns      = (/0, 1, 0, 1, 0, 1 /) ; 0, 0, 0, 0/)
res@xyLineColors        = (/"blue","blue","red","red", "green", "green"/)
                      ; "darkgreen","purple","blue","chartreuse","orange","cyan"/)          ; change line color
res@xyExplicitLegendLabels = case_name

; res1@cnFillMode       = "RasterFill"            ; Raster Mode
res1@cnFillOn           = True                            ;��ɫ
res1@cnLinesOn          = False                       ; �رյ�ֵ������
res1@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"


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

x             = ispan(1,time_tot,1)

; ______________________ plot top-mid-bottom soil moisture _____________________
plots            = new(4,graphic)
res@trYMinF      = 0
res@trYMaxF      = 0.45

print("plot SM64cm")
res@gsnLeftString= "Soil moisture in the top 64cm (m3/m3)"
SM64cm_arry      = new ( (/4,time_tot/), typeof(SM64cm))
SM64cm_rc_ctl    = regline(x, SM64cm(0,:))          ; slope
SM64cm_rc_sen    = regline(x, SM64cm(1,:))          ; slope
SM64cm_arry(0,:) = SM64cm(0,:)
SM64cm_arry(1,:) = SM64cm_rc_ctl*x + SM64cm_rc_ctl@yintercept
SM64cm_arry(2,:) = SM64cm(1,:)
SM64cm_arry(3,:) = SM64cm_rc_sen*x + SM64cm_rc_sen@yintercept
plots(0)         = gsn_csm_xy(wks,x,SM64cm_arry,res)

print("plot 4th SM")
res@gsnLeftString  = "Soil moisture in 64-173cm (m3/m3)"
SM4th_arry      = new ( (/4,time_tot/), float)
SM4th_rc_ctl    = regline(x, SoilMoist(0,4,:))          ; slope
SM4th_rc_sen    = regline(x, SoilMoist(1,4,:))          ; slope
SM4th_arry(0,:) = SoilMoist(0,4,:)
SM4th_arry(1,:) = SM4th_rc_ctl*x + SM4th_rc_ctl@yintercept
SM4th_arry(2,:) = SoilMoist(1,4,:)
SM4th_arry(3,:) = SM4th_rc_sen*x + SM4th_rc_sen@yintercept
plots(1)        = gsn_csm_xy(wks,x,SM4th_arry,res)

print("plot 5th SM")
res@gsnLeftString  = "Soil moisture in 173-460cm (m3/m3)"
SM5th_arry      = new ( (/4,time_tot/), float)
SM5th_rc_ctl    = regline(x, SoilMoist(0,5,:))          ; slope
SM5th_rc_sen    = regline(x, SoilMoist(1,5,:))          ; slope
SM5th_arry(0,:) = SoilMoist(0,5,:)
SM5th_arry(1,:) = SM5th_rc_ctl*x + SM5th_rc_ctl@yintercept
SM5th_arry(2,:) = SoilMoist(1,5,:)
SM5th_arry(3,:) = SM5th_rc_sen*x + SM5th_rc_sen@yintercept
plots(2)        = gsn_csm_xy(wks,x,SM5th_arry,res)

res@gsnLeftString = "Soil Moisture Difference"
delete(res@trYMinF)
delete(res@trYMaxF)
diff      = new ( (/6,time_tot/), float)
diff(0,:) = SM64cm_arry(2,:) - SM64cm_arry(0,:)
diff(1,:) = SM64cm_arry(3,:) - SM64cm_arry(1,:)
diff(2,:) = SM4th_arry(2,:) - SM4th_arry(0,:)
diff(3,:) = SM4th_arry(3,:) - SM4th_arry(1,:)
diff(4,:) = SM5th_arry(2,:) - SM5th_arry(0,:)
diff(5,:) = SM5th_arry(3,:) - SM5th_arry(1,:)
plots(3)  = gsn_csm_xy(wks,x,diff,res)

; ; plot discharge & recharge
; res@gsnLeftString  = "Qrec: discharge = 1; recharge = -1"
; Qrec(2,:)          = 0.
; plot               = gsn_csm_xy(wks,x,Qrec,res)
; overlay(plot,plots(3))

gsn_panel(wks,(/plots/),(/4,1/),pres)
delete(plots)

; _______________________________ plot soil profile _____________________________
plots       = new(6,graphic)
SoilMoist!0 = "case"
SoilMoist!1 = "depth"
SoilMoist!2 = "time"
SoilMoist&depth = Depth
SoilMoist&time  = x
res@gsnMaximize         = False
pres@gsnPanelLabelBar   = True
res1@tmYLMode           = "Explicit"
res1@tmBorderThicknessF = 3.0
res1@tmYLValues         = ispan(1,7,1)
res1@tiXAxisString      = "Depth (cm)"
res1@tmYLLabels         = tostring(toint(Depth(0:6:-1)))

res2 = res
res2@gsnLeftString       = "Soil Moisture Profile"
res2@vpHeightF           = 0.36
res2@vpWidthF            = 0.2
res2@gsnPaperOrientation = "portrait";"portrait";"landscape" ;"auto"
res2@tmYLMode            = "Explicit"
delete(res2@tmXBValues)
delete(res2@tmXBLabels)
res2@tmYLValues         = ispan(1,7,1)
res2@tiYAxisString      = "Depth (cm)"
res2@tmYLLabels         = tostring(toint(Depth(0:6:-1)))
res2@tmBorderThicknessF = 2.0
;res2@trXMinF            = 0
;res2@trXMaxF            = 0.40
delete(res2@tiXAxisString)

res1@cnLevels      = fspan(0.0,0.4,41)
res1@cnFillPalette = "GMT_drywet" ;"WhiteBlueGreenYellowRed"
res1@gsnLeftString = case_name(0)
plots(0)           = gsn_csm_contour(wks, SoilMoist(0,0:6:-1,:),res1)
plots(1)           = gsn_csm_xy(wks,dim_avg_n_Wrap(SoilMoist(0,0:6:-1,:),1),ispan(1,7,1),res2)

res1@gsnLeftString = case_name(1)
plots(2)           = gsn_csm_contour(wks, SoilMoist(1,0:6:-1,:),res1)
plots(3)           = gsn_csm_xy(wks,dim_avg_n_Wrap(SoilMoist(1,0:6:-1,:),1),ispan(1,7,1),res2)
delete(res1@cnLevels)

res1@cnLevels      = fspan(-0.1,0.1,21)
res1@cnFillPalette = "MPL_BrBG" ;"ViBlGrWhYeOrRe"
res1@gsnLeftString =  case_name(1) +" - "+  case_name(0)
plots(4)           = gsn_csm_contour(wks, SoilMoist(1,0:6:-1,:) - SoilMoist(0,0:6:-1,:),res1)
;res2@trXMinF       = -0.1
;res2@trXMaxF       = 0.1
plots(5)           = gsn_csm_xy(wks,dim_avg_n_Wrap(SoilMoist(:,0:6:-1,:),2),ispan(1,7,1),res2)
gsn_panel(wks,(/plots/),(/3,2/),pres)

delete(plots)

; _____________________________ VegT TVeg Fwsoil  ______________________________
res@gsnMaximize         = True
pres@gsnPanelLabelBar   = False
plots = new(6,graphic)

; plot VegT
print("plot vegetation temperature")
res@gsnLeftString  = "vegetation temperature (C)"
res@trYMinF    = 0.
res@trYMaxF    = 45.
VegT_arry      = new ( (/4,time_tot/), typeof(VegT))
VegT_rc_ctl    = regline(x, VegT(0,:))          ; slope
VegT_rc_sen    = regline(x, VegT(1,:))          ; slope
VegT_arry(0,:) = VegT(0,:)
VegT_arry(1,:) = VegT_rc_ctl*x + VegT_rc_ctl@yintercept
VegT_arry(2,:) = VegT(1,:)
VegT_arry(3,:) = VegT_rc_sen*x + VegT_rc_sen@yintercept
plots(0) = gsn_csm_xy(wks,x,VegT_arry,res)

res@gsnLeftString  = "vegetation temperature diff (C)"
res@trYMinF        = -1.
res@trYMaxF        = 1.
plots(1) = gsn_csm_xy(wks,x,VegT(1,:) - VegT(0,:),res)


; plot TVeg
print("plot vegetation temperature")
res@gsnLeftString  = "Transpiration (mm/d)"
res@trYMinF        = 0.
res@trYMaxF        = 3.
TVeg_arry      = new ( (/4,time_tot/), typeof(TVeg))
TVeg_rc_ctl    = regline(x, TVeg(0,:))          ; slope
TVeg_rc_sen    = regline(x, TVeg(1,:))          ; slope
TVeg_arry(0,:) = TVeg(0,:)
TVeg_arry(1,:) = TVeg_rc_ctl*x + TVeg_rc_ctl@yintercept
TVeg_arry(2,:) = TVeg(1,:)
TVeg_arry(3,:) = TVeg_rc_sen*x + TVeg_rc_sen@yintercept
plots(2) = gsn_csm_xy(wks,x,TVeg_arry,res)

res@gsnLeftString  = "Transpiration diff (mm/d)"
res@trYMinF        = -0.5
res@trYMaxF        = 0.5
plots(3) = gsn_csm_xy(wks,x,TVeg(1,:) - TVeg(0,:),res)

; plot Fwsoil
print("plot Fwsoil")
res@gsnLeftString  = "Fwsoil (-)"
res@trYMinF        = 0.
res@trYMaxF        = 1.
Fwsoil_arry      = new ( (/4,time_tot/), typeof(Fwsoil))
Fwsoil_rc_ctl    = regline(x, Fwsoil(0,:))          ; slope
Fwsoil_rc_sen    = regline(x, Fwsoil(1,:))          ; slope
Fwsoil_arry(0,:) = Fwsoil(0,:)
Fwsoil_arry(1,:) = Fwsoil_rc_ctl*x + Fwsoil_rc_ctl@yintercept
Fwsoil_arry(2,:) = Fwsoil(1,:)
Fwsoil_arry(3,:) = Fwsoil_rc_sen*x + Fwsoil_rc_sen@yintercept
plots(4) = gsn_csm_xy(wks,x,Fwsoil_arry,res)

res@gsnLeftString  = "Fwsoil diff (-)"
res@trYMinF        = -0.2
res@trYMaxF        = 0.2
plots(5) = gsn_csm_xy(wks,x,Fwsoil(1,:) - Fwsoil(0,:),res)
gsn_panel(wks,(/plots/),(/3,2/),pres)
delete(res@trYMinF)
delete(res@trYMaxF)
delete(plots)
;
; ; _________________________________ EVAPFBL ____________________________________
; res@trYMinF = 0.
; res@trYMaxF = 1.2
;
; plots = new(4,graphic)
; print("plot EVAP64cm")
; res@gsnLeftString  = "Transpiration taken from the top 64cm (mm/d)"
; EVAP64cm_arry      = new ( (/4,time_tot/), typeof(EVAP64cm))
; EVAP64cm_rc_ctl    = regline(x, EVAP64cm(0,:))          ; slope
; EVAP64cm_rc_sen    = regline(x, EVAP64cm(1,:))          ; slope
; EVAP64cm_arry(0,:) = EVAP64cm(0,:)
; EVAP64cm_arry(1,:) = EVAP64cm_rc_ctl*x + EVAP64cm_rc_ctl@yintercept
; EVAP64cm_arry(2,:) = EVAP64cm(1,:)
; EVAP64cm_arry(3,:) = EVAP64cm_rc_sen*x + EVAP64cm_rc_sen@yintercept
; plots(0) = gsn_csm_xy(wks,x,EVAP64cm_arry,res)
;
; print("plot 4th EVAPFBL")
; res@gsnLeftString  = "Transpiration taken from 64-173cm (mm/d)"
; EVAP4th_arry      = new ( (/4,time_tot/), float)
; EVAP4th_rc_ctl    = regline(x, EVAPFBL(0,4,:))          ; slope
; EVAP4th_rc_sen    = regline(x, EVAPFBL(1,4,:))          ; slope
; EVAP4th_arry(0,:) = EVAPFBL(0,4,:)
; EVAP4th_arry(1,:) = EVAP4th_rc_ctl*x + EVAP4th_rc_ctl@yintercept
; EVAP4th_arry(2,:) = EVAPFBL(1,4,:)
; EVAP4th_arry(3,:) = EVAP4th_rc_sen*x + EVAP4th_rc_sen@yintercept
; plots(1) = gsn_csm_xy(wks,x,EVAP4th_arry,res)
;
; print("plot 5th EVAPBL")
; res@gsnLeftString  = "Transpiration taken from 173-460cm (mm/d)"
; EVAP5th_arry      = new ( (/4,time_tot/), float)
; EVAP5th_rc_ctl    = regline(x, EVAPFBL(0,5,:))          ; slope
; EVAP5th_rc_sen    = regline(x, EVAPFBL(1,5,:))          ; slope
; EVAP5th_arry(0,:) = EVAPFBL(0,5,:)
; EVAP5th_arry(1,:) = EVAP5th_rc_ctl*x + EVAP5th_rc_ctl@yintercept
; EVAP5th_arry(2,:) = EVAPFBL(1,5,:)
; EVAP5th_arry(3,:) = EVAP5th_rc_sen*x + EVAP5th_rc_sen@yintercept
; plots(2) = gsn_csm_xy(wks,x,EVAP5th_arry,res)
;
;
; print("plot EVAPFBL_diff")
; res@gsnLeftString  = "EVAPFBL_diff (SEN-CTL) (mm/d)"
; res@trYMinF        = -0.5
; res@trYMaxF        = 0.5
; EVAPFBL_diff       = new((/3,time_tot/), float)
; delete(res@xyLineThicknesses)
; delete(res@xyDashPatterns)
; delete(res@xyLineColors)
; res@xyLineThicknesses   = (/  1.5, 1.5, 1.5/)
; res@xyDashPatterns      = (/0, 0, 0/)
; res@xyLineColors        = (/"blue","red","green"/)
; EVAPFBL_diff(0,:) = EVAP64cm(1,:) - EVAP64cm(0,:)
; EVAPFBL_diff(1,:) = EVAPFBL(1,4,:) - EVAPFBL(0,4,:)
; EVAPFBL_diff(2,:) = EVAPFBL(1,5,:) - EVAPFBL(0,5,:)
; plots(3) = gsn_csm_xy(wks,x,EVAPFBL_diff,res)
; gsn_panel(wks,(/plots/),(/4,1/),pres)

end
