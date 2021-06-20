
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
;*******************************************

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 500000000000 ;
end setvalues
;_____________________________________

; =========================== Notes for time step data =========================
; model time step : 0, 1, 2, 3 , 4 , 5,  6,  7
; real world time : 0, 3, 6, 9, 12, 15, 18, 21
; ==============================================================================

undef("add_x_errorbar")
procedure add_x_errorbar(wks, plot, x[*]:numeric,y[*]:numeric, xerror[*]:numeric, res)
local dimz, xpoints, ypoints, gsres, i
begin
  gsres = res
  dimz  = dimsizes(y)

  xpoints = (/x-xerror, x+xerror/)
  ypoints = (/ y, y/)

  do i=0, dimz(0)-1
     plot@$tostring("x_errorbar"+i)$ = gsn_add_polyline(wks, plot, xpoints(:,i), ypoints(:,i), gsres )
  end do
end
;---
undef("add_y_errorbar")
procedure add_y_errorbar(wks, plot, x[*]:numeric,y[*]:numeric, yerror[*]:numeric, res)
local dimz, xpoints, ypoints, gsres, i
begin
  gsres = res
  dimz  = dimsizes(y)

  xpoints = (/x, x/)
  ypoints = (/ y-yerror, y+yerror/)
  do i=0, dimz(0)-1
     plot@$tostring("y_errorbar"+i)$ = gsn_add_polyline(wks, plot, xpoints(:,i), ypoints(:,i), gsres )
  end do
end



begin

; =========================== OPTIONS ===========================
PFT        = "tree" ; "tree"
scale      = "SE-AU"
year_s     = 2019
year_e     = 2019

vpd_fwsoil = True

HW_day     = (/ 15, 25 /)
event_sum  = dimsizes(HW_day)

message    = "LST_subday"

pic_name  = message+"_"+scale+"_"+year_s+"-"+year_e

path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"

case_name = (/"CTL-30x3+30yr","FREE_DRAIN","CTL-30x3+30yr_2019-DeepRoot"/)
case_sum  = dimsizes(case_name)
case_name_brf = (/"GW","FD","DR"/)

; =========================== Pre-define ===========================
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

; sub region 1
A_rgn_1 = -36
B_rgn_1 = -34
C_rgn_1 = 149
D_rgn_1 = 151

; sub region 2
A_rgn_2 = -38
B_rgn_2 = -36
C_rgn_2 = 146.5
D_rgn_2 = 148.5

; ______ load path ______
path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = path_mum+case_name(case_num)+"/outputs-subday"
end do
filename1 = path(0)+"/cable_out_2019_SE_Aus.nc"
f1        = addfile (filename1,"r")
lat       = f1->y({A:B})
lon       = f1->x({C:D})
latitude  = f1->latitude({A:B},{C:D})
longitude = f1->longitude({A:B},{C:D})
lat_sum   = dimsizes(lat)
lon_sum   = dimsizes(lon)
default_val = f1->Qle@_FillValue
delete(f1)
delete(filename1)

order_name  = (/"(a)","(b)","(c)","(d)",  "(e)","(f)","(g)","(h)",\
                "(i)","(j)","(k)","(l)",  "(m)","(n)","(o)","(p)",\
                "(q)","(r)","(s)","(t)"/)

; ============================= testing plot ===================================
wks1 = gsn_open_wks("pdf","Check_LST")
gsn_define_colormap(wks1,"rainbow")

res9                         = True
res9@gsnDraw                 = True
res9@gsnFrame                = True
res9@gsnMaximize             = False
res9@cnFillMode              = "RasterFill"            ; Raster Mode
res9@cnFillOn                = True
res9@cnLinesOn               = False

; ================================== VegT ======================================
; Note that:
;    Even VegT & Tair are set to hold 31 days,
;    actually only 15th and 25th Jan are given to VegT & Tair, and other days are 0

VegT     = new((/case_sum,2,31,lat_sum,lon_sum/),float,default_val) ; case x Terra/Aqua x days x lat x lat
Tair     = new((/2,31,lat_sum,lon_sum/),float,default_val)          ; Terra/Aqua x days x lat x lon
VegT     = 0.
Tair     = 0.

do case_num = 0, case_sum -1
  print("case_num = " +case_num)
  do year = year_s,year_e
    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
    f1 = addfile (filename1,"r")

    VegT_tmp   = f1->VegT(0:60*8,{A:B},{C:D})
    do event_num = 0,event_sum-1
      ; 0 am
      ts = (HW_day(event_num)-1)*8
      ; to compare with Terra - 10am
      VegT(case_num, 0, HW_day(event_num)-1, :, :) = VegT_tmp(ts+3,:,:) + (VegT_tmp(ts+4,:,:) - VegT_tmp(ts+3,:,:))/3
      ; to compare with Aqua  - 2pm
      VegT(case_num, 1, HW_day(event_num)-1, :, :) = VegT_tmp(ts+4,:,:) + (VegT_tmp(ts+5,:,:) - VegT_tmp(ts+4,:,:))*(2/3)
    end do

    if case_num .eq. 0 then
      Tair_tmp   = f1->Tair(0:60*8,{A:B},{C:D})
      do event_num = 0,event_sum-1
        ; 0 am
        ts = (HW_day(event_num)-1)*8
        ; to compare with Terra - 10am
        Tair(0,HW_day(event_num)-1,:,:) =  Tair_tmp(ts+3,:,:) + (Tair_tmp(ts+4,:,:) - Tair_tmp(ts+3,:,:))/3
        ; to compare with Aqua  - 2pm
        Tair(1,HW_day(event_num)-1,:,:) =  Tair_tmp(ts+4,:,:) + (Tair_tmp(ts+5,:,:) - Tair_tmp(ts+4,:,:))*(2/3)
      end do
    end if
  end do ; year
end do ; case_num

VegT        = VegT - 273.15
Tair        = Tair - 273.15

VegT!0      = "case"
VegT!1      = "time_step"
VegT!2      = "time"
VegT!3      = "lat"
VegT!4      = "lon"
VegT&lat    = lat
VegT&lon    = lon

plot = gsn_csm_contour(wks1,VegT(0,0,14,:,:),res9)
plot = gsn_csm_contour(wks1,VegT(0,1,14,:,:),res9)

; ================================== LST =======================================
path_LST = "/g/data/w35/mm3972/data/MODIS_LST/";
Terra    = "MOD11A1.006_1km_aid0001.nc"
Aqua     = "MYD11A1.006_1km_aid0001.nc"

filename2 = path_LST+Terra
f2  = addfile (filename2,"r")

filename3 = path_LST+Aqua
f3  = addfile (filename3,"r")

; ______________________ Terra _______________________
; Terra passes over Aus ~ 10 am
LST_Terra_QC = f2->QC_Day(:,{A:B},{C:D})
; LST_Terra_QC = f2->QC_Night(:,{A:B},{C:D})
printVarSummary(LST_Terra_QC)

LST_Terra = tofloat(f2->LST_Day_1km(:,{A:B},{C:D}))
; LST_Terra = tofloat(f2->LST_Night_1km(:,{A:B},{C:D}))
LST_Terra = where((LST_Terra .ne. 0) .and. (LST_Terra_QC .eq. 0) ,\
                 LST_Terra*0.02-273.15,-9999.)
LST_Terra@_FillValue = -9999.

LST_Terra_regrid = linint2_Wrap(doubletofloat(f2->lon({C:D})),doubletofloat(f2->lat({A:B})),LST_Terra,False,lon,lat,0)
LST_Terra_regrid!0 = "time"
LST_Terra_regrid!1 = "lat"
LST_Terra_regrid!2 = "lon"
LST_Terra_regrid&lat = lat
LST_Terra_regrid&lon = lon
LST_Terra_regrid@_FillValue = -9999.
; printVarSummary(LST_Terra_regrid)


res9@tiYAxisString = "LST_Terra_regrid(14,:,:)"
plot = gsn_csm_contour(wks1,LST_Terra_regrid(14,:,:),res9)
res9@tiYAxisString = "LST_Terra_regrid(24,:,:)"
plot = gsn_csm_contour(wks1,LST_Terra_regrid(24,:,:),res9)

; ______________________ Aqua _______________________
; Aqua passes over Aus ~ 2 pm
LST_Aqua_QC = f3->QC_Day(:,{A:B},{C:D})
; LST_Aqua_QC = f3->QC_Night(:,{A:B},{C:D})
printVarSummary(LST_Aqua_QC)

LST_Aqua = tofloat(f3->LST_Day_1km(:,{A:B},{C:D}))
; LST_Aqua = tofloat(f3->LST_Night_1km(:,{A:B},{C:D}))
LST_Aqua = where((LST_Aqua .ne. 0) .and. (LST_Aqua_QC .eq. 0) ,\
                LST_Aqua*0.02-273.15,-9999.)
LST_Aqua@_FillValue = -9999.

LST_Aqua_regrid = linint2_Wrap(doubletofloat(f3->lon({C:D})),doubletofloat(f3->lat({A:B})),LST_Aqua,False,lon,lat,0)
LST_Aqua_regrid!0 = "time"
LST_Aqua_regrid!1 = "lat"
LST_Aqua_regrid!2 = "lon"
LST_Aqua_regrid&lat = lat
LST_Aqua_regrid&lon = lon
LST_Aqua_regrid@_FillValue = -9999.
; printVarSummary(LST_Aqua_regrid)

res9@tiYAxisString = "LST_Aqua_regrid(14,:,:)"
plot = gsn_csm_contour(wks1,LST_Aqua_regrid(14,:,:),res9)
res9@tiYAxisString = "LST_Aqua_regrid(24,:,:)"
plot = gsn_csm_contour(wks1,LST_Aqua_regrid(24,:,:),res9)

; ================================ Calc Diff ===================================
VegT(0,:,:,:,:) = VegT(0,:,:,:,:) - Tair(:,:,:,:)
VegT(1,:,:,:,:) = VegT(1,:,:,:,:) - Tair(:,:,:,:)
VegT(2,:,:,:,:) = VegT(2,:,:,:,:) - Tair(:,:,:,:)
LST_Terra_regrid= LST_Terra_regrid - Tair(0,:,:,:)
LST_Aqua_regrid = LST_Aqua_regrid  - Tair(1,:,:,:)


; =============================== Diurnal Cycle ================================

; ---------------- Mask ----------------
; mask trees
filename1 = path(0)+"/cable_out_2019_SE_Aus.nc"
f1        = addfile (filename1,"r")

PFT_mask        = new((/3,lat_sum,lon_sum/),float,f1->iveg@_FillValue) ; tree/region1/region2 x lat x lon

; tree region
PFT_mask(0,:,:) = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, -1)
; tree region1
PFT_mask(1,:,:) = where(  \
                        (f1->latitude({A:B},{C:D}) .ge. A_rgn_1)  .and. (f1->latitude({A:B},{C:D}) .le. B_rgn_1)   \
                  .and. (f1->longitude({A:B},{C:D}) .ge. C_rgn_1) .and. (f1->longitude({A:B},{C:D}) .le. D_rgn_1), \
                  PFT_mask(0,:,:),-1)
; tree region2
PFT_mask(2,:,:) = where(  \
                        (f1->latitude({A:B},{C:D}) .ge. A_rgn_2)  .and. (f1->latitude({A:B},{C:D}) .le. B_rgn_2)   \
                  .and. (f1->longitude({A:B},{C:D}) .ge. C_rgn_2) .and. (f1->longitude({A:B},{C:D}) .le. D_rgn_2),  \
                  PFT_mask(0,:,:), -1)

PFT_mask_multi = new((/9,3,lat_sum,lon_sum/),float,PFT_mask@_FillValue) ; timestep x tree/region1/region2 x lat x lon
do z = 0, 8
  PFT_mask_multi(z,:,:,:) = PFT_mask
end do

res9@tiYAxisString = "PFT_mask(1,:,:)"
plot = gsn_csm_contour(wks1,PFT_mask(1,:,:),res9)
res9@tiYAxisString = "PFT_mask(2,:,:)"
plot = gsn_csm_contour(wks1,PFT_mask(2,:,:),res9)

delete(filename1)
delete(f1)

; HW mask
filename_hw = "/g/data/w35/mm3972/scripts/ehfheatwaves/nc_file/SE_AUS/HW_Event_Indicator_1970-2019.nc" ; summer HW only
fhw         = addfile(filename_hw,"r")

HW = new((/31,lat_sum,lon_sum/),float)  ; timestep x lat x lon
day_start = 0
do year = 1970, 2019 - 1
  if ( mod(year,4) .eq. 0) then
     day_start = day_start + 366
  else
     day_start = day_start + 365
  end if
end do
HW = fhw->event(day_start:day_start+31-1,:,:)*1.
delete(filename_hw)
delete(fhw)

HW_mask_multi = new((/2,9,lat_sum,lon_sum/),float) ; hw x timestep x lat x lon
do z = 0, 8
  ; 15 Jan 2019
  HW_mask_multi(0,z,:,:) = HW(14,:,:)
  ; 25 Jan 2019
  HW_mask_multi(1,z,:,:) = HW(24,:,:)
end do
; delete(HW)

res9@tiYAxisString = "HW_mask_multi(0,0,:,:)"
plot = gsn_csm_contour(wks1,HW_mask_multi(0,0,:,:),res9)
res9@tiYAxisString = "HW_mask_multi(1,4,:,:)"
plot = gsn_csm_contour(wks1,HW_mask_multi(1,4,:,:),res9)

; LST
LST_mask = new((/31,9,lat_sum,lon_sum/),float,-9999.) ; days x timestep x lat x lon
LST_mask = 0.
do z = 0, 8
  LST_mask(:,z,:,:) = where(LST_Terra_regrid .ne. -9999. .and. LST_Aqua_regrid .ne. -9999., 1., 2.)
end do
res9@tiYAxisString = "LST_mask(14,0,:,:)"
plot = gsn_csm_contour(wks1,LST_mask(14,0,:,:),res9)
res9@tiYAxisString = "LST_mask(24,0,:,:)"
plot = gsn_csm_contour(wks1,LST_mask(24,0,:,:),res9)

; ---------------- Diurnal Cycle ----------------
deltaT   = new((/case_sum,2,3,9/), float, default_val) ; case x hw x tree/tree region1/tree region2 x time steps
deltaT_u = new((/case_sum,2,3,9/), float, default_val)
deltaT   = 0.
deltaT_u = 0.

do case_num = 0, case_sum-1
  print("case_num = " +case_num)
  filename1 = path(case_num)+"/cable_out_2019_SE_Aus.nc"
  f1 = addfile (filename1,"r")

  VegT_tmp   = f1->VegT(0:60*8,{A:B},{C:D})

  do hw = 0, dimsizes(HW_day) - 1
    ts = (HW_day(hw)-1)*8
    do region = 0,2 ; tree/tree region1/tree region2
      ; mean
      deltaT(case_num,hw,region,:) = dim_avg_n_Wrap(where(\
                                              LST_mask(HW_day(hw)-1,:,:,:) .eq. 1. \
                                              .and. PFT_mask_multi(:,region,:,:) .eq. 1 \
                                              .and. HW_mask_multi(hw,:,:,:) .eq. 1, \
                                              VegT_tmp(ts:ts+8,:,:) - Tair_tmp(ts:ts+8,:,:), default_val) ,(/1,2/))
      ; Uncertainty
      deltaT_u(case_num,hw,region,:) = dim_stddev_n_Wrap(where(\
                                              LST_mask(HW_day(hw)-1,:,:,:) .eq. 1. \
                                              .and. PFT_mask_multi(:,region,:,:) .eq. 1 \
                                              .and. HW_mask_multi(hw,:,:,:) .eq. 1, \
                                              VegT_tmp(ts:ts+8,:,:) - Tair_tmp(ts:ts+8,:,:), default_val) ,(/1,2/))
    end do
  end do
  delete(filename1)
  delete(f1)
end do ; case_num

; uncertainty
deltaT_uncert = new((/case_sum,2,3,2,9/), float, default_val)  ; case x hw x tree/tree region1/tree region2 x up/low boundary x time steps
deltaT_uncert(:,:,:,0,:) = deltaT-deltaT_u
deltaT_uncert(:,:,:,1,:) = deltaT+deltaT_u

; ----------------- LST dots ------------------
LST_tree   = new((/2,3,2/),float,-9999.) ; hw x regions x 2 time steps
LST_tree_u = new((/2,3,2/),float,-9999.) ; hw x regions x 2 time steps

do region = 0,2
  ; mean
  LST_tree(0,region,0) = avg(where( PFT_mask(region,:,:) .eq. 1 .and. HW(14,:,:) .eq. 1\
                             .and. LST_Aqua_regrid(14,:,:) .ne. -9999., LST_Terra_regrid(14,:,:), -9999.) )
  LST_tree(1,region,0) = avg(where( PFT_mask(region,:,:) .eq. 1 .and. HW(24,:,:) .eq. 1\
                             .and. LST_Aqua_regrid(24,:,:) .ne. -9999., LST_Terra_regrid(24,:,:), -9999.) )
  LST_tree(0,region,1) = avg(where( PFT_mask(region,:,:) .eq. 1 .and. HW(14,:,:) .eq. 1\
                             .and. LST_Terra_regrid(14,:,:) .ne. -9999., LST_Aqua_regrid(14,:,:), -9999.) )
  LST_tree(1,region,1) = avg(where( PFT_mask(region,:,:) .eq. 1 .and. HW(24,:,:) .eq. 1\
                             .and. LST_Terra_regrid(24,:,:) .ne. -9999., LST_Aqua_regrid(24,:,:), -9999.) )
  ; Uncertainty
  LST_tree_u(0,region,0) = stddev(where( PFT_mask(region,:,:) .eq. 1 .and. HW(14,:,:) .eq. 1\
                             .and. LST_Aqua_regrid(14,:,:) .ne. -9999., LST_Terra_regrid(14,:,:), -9999.))
  LST_tree_u(1,region,0) = stddev(where( PFT_mask(region,:,:) .eq. 1 .and. HW(24,:,:) .eq. 1\
                            .and. LST_Aqua_regrid(24,:,:) .ne. -9999., LST_Terra_regrid(24,:,:), -9999.))
  LST_tree_u(0,region,1) = stddev(where( PFT_mask(region,:,:) .eq. 1 .and. HW(14,:,:) .eq. 1\
                             .and. LST_Terra_regrid(14,:,:) .ne. -9999., LST_Aqua_regrid(14,:,:), -9999.))
  LST_tree_u(1,region,1) = stddev(where( PFT_mask(region,:,:) .eq. 1 .and. HW(24,:,:) .eq. 1\
                            .and. LST_Terra_regrid(24,:,:) .ne. -9999., LST_Aqua_regrid(24,:,:), -9999.))
end do


; ================================= Plotting ===================================
pic = "./plots/Fig6_revision_"+pic_name
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")

res                         = True
res@gsnDraw                 = False  ; Don't draw plots
res@gsnFrame                = False  ; ����ҳ
res@gsnMaximize             = False

res@cnFillOn                = True
res@cnMissingValFillColor   = "gray80"
res@cnFillBackgroundColor   = "gray80"

res@vpHeightF 	            = 0.38               ; change aspect ratio of plot
res@vpWidthF 	              = 0.38
res@tiYAxisFontThicknessF   = 3.
res@tmXBMajorThicknessF     = 3.
res@tmXBMinorThicknessF     = 3.
res@tmYLMajorThicknessF     = 3.
res@tmYLMinorThicknessF     = 3.
res@tmYLLabelFontThicknessF = 3.
res@tmBorderThicknessF      = 3.
; res@txFont                  = 25
res@txFontHeightF           = 0.015            ; change title font heights
res@gsnPaperOrientation     = "landscape";"portrait";"landscape" ;"auto"

res@cnFillMode             = "RasterFill"            ; Raster Mode
res@cnFillOn               = True                            ;��ɫ
res@cnLinesOn              = False
res@cnLineLabelsOn         = False
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode   = "ExplicitLevels"

res@lbLabelBarOn           = True
res@lbBoxLinesOn           = True                       ;�ر�lbar box ����
res@lbOrientation          = "vertical";"vertical" ;"horizontal"
res@lbTitleFont            = 25
res@lbLabelFont            = 25
res@lbTitleFontHeightF     = 0.02
res@lbLabelFontHeightF     = 0.015

res@tmXBLabelFont           = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont           = 25
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF    = 0.015
res@tmYLLabelFontHeightF    = 0.015
; res@mpDataBaseVersion       = "HighRes"

res@lbLabelBarOn            = True
res@lbTitlePosition         = "Right"
res@lbTitleAngleF           = 270
res@lbTitleDirection        = "Across"                          ; letter angle
res@tmYLMode                = "Explicit"
res@tmXBMode                = "Explicit"
res@txString                = ""
res@gsnLeftString           = ""
res@gsnRightString          = ""
if scale .eq. "Global" then
  res@trYMinF         = -60
  res@trYMaxF         = 80
  res@mpMinLatF       = -60
  res@mpMaxLatF       = 80
  res@mpCenterLonF    = 160
  res@tmYLMode        = "Explicit"
  res@tmYLValues      =(/60,30,0,-30,-60/)
  res@tmYLLabels      =(/"60~S~o~N~N","30~S~o~N~N","0~S~o~N~","30~S~o~N~S","60~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      =(/240,300,0,60,120,180/)
  res@tmXBLabels      =(/"120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/)
else if scale .eq. "SE-AU" then
  res@gsnAddCyclic    = False
  res@trYMinF         = -40.
  res@trYMaxF         = -28.
  res@mpMinLatF       = -40.
  res@mpMaxLatF       = -28.
  res@trXMinF         = 140.
  res@trXMaxF         = 154.
  res@mpMinLonF       = 140.
  res@mpMaxLonF       = 154.
  res@tmYLMode        = "Explicit"
  res@tmYLValues      = (/-40,-35,-30,-25/)
  res@tmYLLabels      = (/"40~S~o~N~S","35~S~o~N~S","30~S~o~N~S","25~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/140,145,150/)
  res@tmXBLabels      = (/"140~S~o~N~E","145~S~o~N~E","150~S~o~N~E"/)
else if scale .eq. "AU" then
  res@gsnAddCyclic    = False
  res@trYMinF         = -44
  res@trYMaxF         = -10
  res@mpMinLatF       = -44
  res@mpMaxLatF       = -10
  res@trXMinF         = 112
  res@trXMaxF         = 154
  res@mpMinLonF       = 112
  res@mpMaxLonF       = 154
  res@tmYLMode        = "Explicit"
  res@tmYLValues      = (/-40,-35,-30,-25,-20,-15,-10/)
  res@tmYLLabels      = (/"40~S~o~N~S","35~S~o~N~S","30~S~o~N~S","25~S~o~N~S","20~S~o~N~S","15~S~o~N~S","10~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/115,120,125,130,135,140,145,150/)
  res@tmXBLabels      = (/"115~S~o~N~E","120~S~o~N~E","125~S~o~N~E","130~S~o~N~E","135~S~o~N~E","140~S~o~N~E","145~S~o~N~E","150~S~o~N~E"/)
end if
end if
end if

; _____________ Colorbar ____________
color_tmp = namedcolor2rgba("brown4");("LightSalmon4")
color_tmp(:,:2) = color_tmp(:,:2) * 255
print(color_tmp) ; -> (/139,35,35,1/)


; Red
cmap1  = (/ (/255.,255,255,1/), (/255,255,229,1/), (/255,247,188,1/), (/254,227,145,1/), \
            (/254,196,79,1/),  (/254,153,41,1/),  (/236,112,20,1/),  (/204,76,2,1/),    \
            (/153,52,4,1/),    (/139,35,35,1/), (/102,37,6,1/)    /)
; Blue
cmap2  = (/ (/255.,255,255,1/), (/247,252,240,1/), (/224,243,219,1/), (/204,235,197,1/),  \
            (/168,221,181,1/), (/123,204,196,1/),  (/78,179,211,1/),  (/43,140,190,1/),  \
            (/8,104,172,1/),    (/  8,64,129,1/),    (/2,56,88,1/) /)

; Blue to Red
cmap3  = (/ (/8.,29,88,1/),    (/  8,64,129,1/),  (/8,104,172,1/),  \                   \
            (/43,140,190,1/),  (/78,179,211,1/),  (/123,204,196,1/),(/168,221,181,1/), \
            (/204,235,197,1/), (/224,243,219,1/), (/255,255,255,1/),\
            (/255,255,255,1/), (/255,247,188,1/), (/254,227,145,1/),\
            (/254,196,79,1/),  (/254,153,41,1/),  (/236,112,20,1/),  (/204,76,2,1/),   \
            (/153,52,4,1/),    (/139,35,35,1/),    (/102,37,6,1/)   /) ;(/128,0,38,1/)   /) ;  (/102,37,6,1/) (/103,0,13,1/) (/2.,56,88,1/),


cmap1(:,0:2) = cmap1(:,0:2)/255.
cmap2(:,0:2) = cmap2(:,0:2)/255.
cmap3(:,0:2) = cmap3(:,0:2)/255.

;**************** ͼ�������� *****************
pres                    = True
pres@wkWidth            = 2500
pres@wkHeight           = 2500
pres@vpHeightF 	        = 0.2               ; change aspect ratio of plot
pres@vpWidthF 	        = 0.95
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False
; pres@gsnPanelMainPosXF  = 0.1
; pres@gsnPanelMainPosYF  = 0.1
pres@cnLevelSelectionMode = "ExplicitLevels"
pres@lbBoxLinesOn       = True
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbTitleFontHeightF = 0.02
pres@lbLabelFontHeightF = 0.02
pres@gsnPanelYWhiteSpacePercent = 0.5
pres@gsnPanelXWhiteSpacePercent = 0.5


; --------- lines in diurnal cycle ---------
res1                    = True
res1@gsnDraw            = True
res1@gsnFrame           = True
res1@gsnMaximize        = False
res1@vpXF               = 0.1               ; change aspect ratio of plot
res1@vpYF 	            = 0.8
res1@vpHeightF          = 0.32               ; change aspect ratio of plot
res1@vpWidthF 	        = 0.38
res1@tiYAxisFontThicknessF   = 3.
res1@tmXBMajorThicknessF     = 3.
res1@tmXBMinorThicknessF     = 3.
res1@tmYLMajorThicknessF     = 3.
res1@tmYLMinorThicknessF     = 3.
res1@tmYLLabelFontThicknessF = 3.
res1@tmBorderThicknessF      = 3.
; res@txFont                  = 25
res1@txFontHeightF           = 0.015            ; change title font heights
res1@gsnPaperOrientation     = "landscape";"portrait";"landscape" ;"auto"


res1@pmLegendWidthF     = 0.10                     ; res1ize legend width
res1@pmLegendHeightF    = 0.10                     ; and height
res1@tiXAxisFontHeightF = 0.02
res1@tiYAxisFontHeightF = 0.02
res1@tmXBLabelFontHeightF = 0.02
res1@tmYLLabelFontHeightF = 0.02

res1@lgPerimOn           = False      ; No legend perimeter.
res1@xyLineThicknesses   = (/2.0, 2.0, 2.0 /)
res1@xyDashPatterns      = (/0, 0, 0/)
res1@xyLineColors        = (/"darkgreen","orange","blue"/)
res1@xyExplicitLegendLabels = case_name_brf
res1@lgItemOrder         = (/ 2,1,0 /)   ; Reorder the legends
res1@lgBoxMinorExtentF      = 0.2
res1@pmLegendParallelPosF   = 0.5
res1@pmLegendOrthogonalPosF = -0.7
res1@pmLegendWidthF         = 0.13
res1@pmLegendHeightF        = 0.11

res1@tmXBMode            = "Explicit"
x                        = ispan(1,9,1)
res1@tmXBValues          = ispan(1,9,1)
res1@tmXBLabels          = tostring(ispan(0,24,3))
res1@trXMinF             = 1
res1@trXMaxF             = 9


; --------- dots in diurnal cycle ---------
res2 = res1
res2@xyMarkLineMode         = "Markers"
res2@xyMarkers              = (/16, 16/)
res2@xyMarkerColors         = (/"red","red"/) ; 3 different colors
delete(res2@xyExplicitLegendLabels)
res2@xyExplicitLegendLabels = (/"Terra","Aqua"/)
res2@pmLegendDisplayMode    = "Never"
res2@lgPerimOn              = False      ; No legend perimeter.
res2@pmLegendParallelPosF   = 0.25
res2@pmLegendOrthogonalPosF = -1.05
res2@lgBoxMinorExtentF      = 0.2
res2@pmLegendWidthF         = 0.20
res2@pmLegendHeightF        = 0.07
res2@xyMarkerSizeF          = 0.008
delete(res2@lgItemOrder)
res2@lgItemOrder            = (/ 1,0 /)   ; Reorder the legend

; -------------- plot order ----------------
txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.03

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

; -------------- start to plot ----------------
plots  = new(10,graphic)
plots2 = new(10,graphic)
plots3 = new(10,graphic)

; -------------- Fig6 a-b LST --------------
res@lbTitleString = "~F5~D~F25~T~B~MOD_2pm~N~ (~S~o~N~C)"
res@cnFillPalette = cmap1 ;"precip3_16lev";"testcmap";"ViBlGrWhYeOrRe";"WhiteBlueGreenYellowRed"
res@cnLevels      = fspan(0,16, 9) ;deltaT_CTL_lvl

res@lbLabelBarOn = False ; "Always"
tmp = LST_Aqua_regrid(14,:,:)
tmp = where(PFT_mask(0,:,:) .eq. 1, tmp, tmp@_FillValue)
plots(0) = gsn_csm_contour_map(wks,tmp,res) ;LST_Aqua_regrid(14,:,:),res)
delete(tmp)

res@lbLabelBarOn = True
tmp = LST_Aqua_regrid(24,:,:)
tmp = where(PFT_mask(0,:,:) .eq. 1, tmp, tmp@_FillValue)
plots(1) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

; -------------- Fig6 c-d deltaT --------------
; notes :
; VegT: case x 2 time steps x days x lat x lon
res@lbTitleString = "~F5~D~F25~T~B~GW_2pm~N~ (~S~o~N~C)"
res@lbLabelBarOn = False ; "Always"
plots(2) = gsn_csm_contour_map(wks,VegT(0,1,14,:,:),res)
; plots(4) = gsn_csm_contour_map(wks,VegT(0,1,15,:,:),res)
res@lbLabelBarOn = True ; "Always"
plots(3) = gsn_csm_contour_map(wks,VegT(0,1,24,:,:),res)
delete(res@cnLevels)
delete(res@cnFillPalette)

; -------------- Fig6 e-f deltaT GW - FD --------------
res@lbTitleString = "~F5~D~F25~T~B~GW_2pm~N~ - ~F5~D~F25~T~B~FD_2pm~N~ (~S~o~N~C)"
res@cnLevels      = (/-5,-4,-3.,-2.,-1,-0.5,0.5,1.,2.,3.,4.,5./)
                   ;fspan(-4,4, 17)
res@cnFillPalette = cmap3;"BlAqGrWh2YeOrReVi22";"BlAqGrWh2YeOrReVi22";"testcmap";"ViBlGrWhYeOrRe";"WhiteBlueGreenYellowRed"
tmp               = VegT(0,1,:,:,:)
tmp               = tmp - VegT(1,1,:,:,:)

res@lbLabelBarOn = False ; "Always"
plots(4) = gsn_csm_contour_map(wks,tmp(14,:,:),res)

res@lbLabelBarOn = True ; "Always"
plots(5) = gsn_csm_contour_map(wks,tmp(24,:,:),res)
delete(tmp)
delete(res@cnLevels)

; -------------- Fig 6 g-h deltaT DR - GW --------------
res@lbTitleString = "~F5~D~F25~T~B~DR_2pm~N~ - ~F5~D~F25~T~B~GW_2pm~N~ (~S~o~N~C)"
res@cnLevels      = (/-5,-4,-3.,-2.,-1,-0.5,0.5,1.,2.,3.,4.,5./)
res@cnFillPalette = cmap3
tmp               = VegT(2,1,:,:,:)
tmp               = tmp - VegT(0,1,:,:,:)

res@lbLabelBarOn = False ; "Always"
plots(6) = gsn_csm_contour_map(wks,tmp(14,:,:),res)

res@lbLabelBarOn = True ; "Always"
plots(7) = gsn_csm_contour_map(wks,tmp(24,:,:),res)
delete(tmp)
delete(res@cnLevels)


; -------------- Fig S6 a-b deltaT GW - LST --------------
res@lbTitleString = "~F5~D~F25~T~B~GW_2pm~N~ - ~F5~D~F25~T~B~MOD_2pm~N~ (~S~o~N~C)"
res@cnLevels      = (/-8,-6,-4,-2.,-1,1.,2.,4.,6.,8./)
res@cnFillPalette = cmap3
tmp               = VegT(0,1,:,:,:)
tmp               = tmp - LST_Aqua_regrid(:,:,:)
tmp1              = tmp(14,:,:)
tmp1              = where(PFT_mask(0,:,:) .eq. 1, tmp1, tmp1@_FillValue)
tmp2              = tmp(24,:,:)
tmp2              = where(PFT_mask(0,:,:) .eq. 1, tmp2, tmp2@_FillValue)
res@lbLabelBarOn = False ; "Always"
plots2(0) = gsn_csm_contour_map(wks,tmp1,res)
res@lbLabelBarOn = True ; "Always"
plots2(1) = gsn_csm_contour_map(wks,tmp2,res)
delete(tmp)
delete(tmp1)
delete(tmp2)

; -------------- Fig S6 c-d deltaT FD - LST --------------
res@lbTitleString = "~F5~D~F25~T~B~FD_2pm~N~ - ~F5~D~F25~T~B~MOD_2pm~N~ (~S~o~N~C)"
res@cnFillPalette = cmap3
tmp               = VegT(1,1,:,:,:)
tmp               = tmp - LST_Aqua_regrid(:,:,:)
tmp1              = tmp(14,:,:)
tmp1              = where(PFT_mask(0,:,:) .eq. 1, tmp1, tmp1@_FillValue)
tmp2              = tmp(24,:,:)
tmp2              = where(PFT_mask(0,:,:) .eq. 1, tmp2, tmp2@_FillValue)
res@lbLabelBarOn = False ; "Always"
plots2(2) = gsn_csm_contour_map(wks,tmp1,res)
res@lbLabelBarOn = True ; "Always"
plots2(3) = gsn_csm_contour_map(wks,tmp2,res)
delete(tmp)
delete(tmp1)
delete(tmp2)
delete(res@cnLevels)

; -------------- Fig S6 e-f deltaT DR - LST --------------
res@lbTitleString = "~F5~D~F25~T~B~DR_2pm~N~ - ~F5~D~F25~T~B~MOD_2pm~N~ (~S~o~N~C)"
res@cnLevels      = (/-8,-6,-4,-2.,-1,1.,2.,4.,6.,8./)
res@cnFillPalette = cmap3
tmp               = VegT(2,1,:,:,:)
tmp               = tmp - LST_Aqua_regrid(:,:,:)
tmp1              = tmp(14,:,:)
tmp1              = where(PFT_mask(0,:,:) .eq. 1, tmp1, tmp1@_FillValue)
tmp2              = tmp(24,:,:)
tmp2              = where(PFT_mask(0,:,:) .eq. 1, tmp2, tmp2@_FillValue)
res@lbLabelBarOn = False ; "Always"
plots2(4) = gsn_csm_contour_map(wks,tmp1,res)
res@lbLabelBarOn = True ; "Always"
plots2(5) = gsn_csm_contour_map(wks,tmp2,res)
delete(tmp)
delete(tmp1)
delete(tmp2)

; -------------- Fig 7 a-b diurnal cycle --------------
; notes:
; deltaT : case x day (0,1,2) x region x time steps
; LST_tree : events x regions x 2 time steps

offset_num = 6
res1@pmLegendDisplayMode = "Always"
res1@trYMinF             = -2.8
res1@trYMaxF             = 13.
y                        = (/4.33, 5.67/) ; 10 am & 2 pm
res1@tiYAxisString       = "~F5~D~F22~T (~S~o~N~C)"
res1@tiXAxisString       = "hour"
res2@tiXAxisString       = "hour"

; ________________ add error bar ________________
gsresy = True
gsresy@gsLineColor      = "red"                       ; add y error bars
gsresy@gsLineThicknessF = 2.0

; region 1
print("region1 15th")
plots3(0) = gsn_csm_xy(wks,x,deltaT(:,0,1,:),res1)
plot1     = gsn_csm_xy(wks,y,LST_tree(0,1,:),res2)
add_y_errorbar(wks, plot1, y, LST_tree(0,1,:), LST_tree_u(0,1,:) ,gsresy)
overlay(plots3(0),plot1)
delete(res1@tiYAxisString)

res1@pmLegendDisplayMode = "Never"

res1@tiXAxisString       = "hour"
print("region2 15th")
plots3(1) = gsn_csm_xy(wks,x,deltaT(:,1,2,:),res1) ; 3+offset_num
plot2     = gsn_csm_xy(wks,y,LST_tree(1,2,:),res2)
add_y_errorbar(wks, plot2, y, LST_tree(1,2,:), LST_tree_u(1,2,:) ,gsresy)
overlay(plots3(1),plot2)

delete(res1@xyLineThicknesses)
delete(res1@xyDashPatterns)
delete(res1@xyLineColors)

print("add uncertainty")
res1@xyLineColor        = -1
res1@gsnXYFillOpacities = 0.1

res1@gsnXYFillColors = "darkgreen"
plot_uncert1  = gsn_csm_xy (wks,x,deltaT_uncert(0,0,1,0:1,:),res1)  ; GW 15th
plot_uncert2  = gsn_csm_xy (wks,x,deltaT_uncert(0,1,2,0:1,:),res1)  ; GW 25th
overlay(plots3(0),plot_uncert1)
overlay(plots3(1),plot_uncert2)

res1@gsnXYFillColors = "orange"
plot_uncert1  = gsn_csm_xy (wks,x,deltaT_uncert(1,0,1,0:1,:),res1)  ; FD 15th
plot_uncert2  = gsn_csm_xy (wks,x,deltaT_uncert(1,1,2,0:1,:),res1)  ; FD 25th
overlay(plots3(0),plot_uncert1)
overlay(plots3(1),plot_uncert2)

res1@gsnXYFillColors = "blue"
plot_uncert1  = gsn_csm_xy (wks,x,deltaT_uncert(2,0,1,0:1,:),res1)  ; DR 15th
plot_uncert2  = gsn_csm_xy (wks,x,deltaT_uncert(2,1,2,0:1,:),res1)  ; DR 25th
overlay(plots3(0),plot_uncert1)
overlay(plots3(1),plot_uncert2)

; ________________ add shading ________________
poly             = new(8,graphic)
gres             = True
gres@gsEdgesOn   = True
gres@gsEdgeColor = "red"
gres@gsLineColor = "red"
gres@gsLineThicknessF = 3.0

xmin = (/C_rgn_1,C_rgn_2/)
xmax = (/D_rgn_1,D_rgn_2/)
ymin = (/A_rgn_1,A_rgn_2/)
ymax = (/B_rgn_1,B_rgn_2/)

xbars = (/ (/xmin(0), xmax(0), xmax(0), xmin(0), xmin(0)/), \
            (/xmin(1), xmax(1), xmax(1), xmin(1), xmin(1)/) /)

ybars = (/ (/ymin(0),ymin(0),ymax(0),ymax(0),ymin(0)/), \
           (/ymin(1),ymin(1),ymax(1),ymax(1),ymin(1)/)  /)

do counter = 0, 7
  if mod(counter,2) .eq. 0 then
    shadow_num = 0
  else
    shadow_num = 1
  end if
  poly(counter) = gsn_add_polyline(wks,plots(counter),xbars(shadow_num,:),ybars(shadow_num,:),gres)
end do

do cnt = 0, 7;9;11
  txid1   = gsn_create_text(wks, order_name(cnt), txres)
  annoid1 = gsn_add_annotation(plots(cnt), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/5,2/),pres) ;(/3,4/),pres)


; _____________ Fig 7 _______________

gsn_panel(wks,(/plots3/),(/5,2/),pres)

; _____________ Fig S6 _______________
poly2 = new(6,graphic)
do cnt = 0, 5
  txid1   = gsn_create_text(wks, order_name(cnt), txres)
  annoid1 = gsn_add_annotation(plots2(cnt), txid1, amres)
end do

do counter = 0, 5
  if mod(counter,2) .eq. 0 then
    shadow_num = 0
  else
    shadow_num = 1
  end if
  poly2(counter) = gsn_add_polyline(wks,plots2(counter),xbars(shadow_num,:),ybars(shadow_num,:),gres)
end do

gsn_panel(wks,(/plots2/),(/5,2/),pres)
delete(plots2)
delete(poly2)


end
