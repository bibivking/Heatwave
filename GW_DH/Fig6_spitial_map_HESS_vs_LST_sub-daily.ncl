
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

begin

; =========================== OPTIONS ===========================
PFT        = "tree" ; "tree"
scale      = "SE-AU"
year_s     = 2019
year_e     = 2019
diff       = "diff"
vpd_fwsoil = True

HW_day_s   = (/ 15, 16, 25 /) ; 15,
HW_day_e   = (/ 15, 16, 25 /) ; 15,
event_sum  = dimsizes(HW_day_s)

day_tot    = 0
do event_num = 0, event_sum-1
  day_tot = day_tot + (HW_day_e(event_num)-HW_day_s(event_num)+1)
end do

message    = "LST_subday_"+diff

pic_name  = "HESS_"+message+"_"+scale+"_"+year_s+"-"+year_e

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
res9@gsnDraw                 = True  ; Don't draw plots
res9@gsnFrame                = True  ; ����ҳ
res9@gsnMaximize             = False
res9@cnFillMode              = "RasterFill"            ; Raster Mode
res9@cnFillOn                = True                            ;��ɫ
res9@cnLinesOn               = False

; ================================== VegT ======================================
VegT     = new((/case_sum,2,31,lat_sum,lon_sum/),float,default_val)
VegT     = 0.
do case_num = 0, case_sum -1
  print("case_num = " +case_num)
  do year = year_s,year_e
    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc" ; outputs/
    f1 = addfile (filename1,"r")

    VegT_tmp   = f1->VegT(0:60*8,{A:B},{C:D})
    do event_num = 0,event_sum-1
      do day = HW_day_s(event_num),HW_day_e(event_num)
        ts = (day-1)*8 ;0 am
        ; to compare with Terra - 10am
        VegT(case_num, 0, day-1, :, :) = VegT_tmp(ts+3,:,:) + (VegT_tmp(ts+4,:,:) - VegT_tmp(ts+3,:,:))/3
        ; to compare with Aqua  - 2pm
        VegT(case_num, 1, day-1, :, :) = VegT_tmp(ts+4,:,:) + (VegT_tmp(ts+5,:,:) - VegT_tmp(ts+4,:,:))*(2/3)
      end do
    end do
  end do ; year
end do ; case_num
; delete(PFT_mask)

VegT        = VegT - 273.15

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
printVarSummary(LST_Terra_QC)

LST_Terra = tofloat(f2->LST_Day_1km(:,{A:B},{C:D}))
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
res9@tiYAxisString = "LST_Terra_regrid(15,:,:)"
plot = gsn_csm_contour(wks1,LST_Terra_regrid(15,:,:),res9)
res9@tiYAxisString = "LST_Terra_regrid(24,:,:)"
plot = gsn_csm_contour(wks1,LST_Terra_regrid(24,:,:),res9)

; ______________________ Aqua _______________________
; Aqua passes over Aus ~ 2 pm
LST_Aqua_QC = f3->QC_Day(:,{A:B},{C:D})
printVarSummary(LST_Aqua_QC)

LST_Aqua = tofloat(f3->LST_Day_1km(:,{A:B},{C:D}))
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
res9@tiYAxisString = "LST_Aqua_regrid(15,:,:)"
plot = gsn_csm_contour(wks1,LST_Aqua_regrid(15,:,:),res9)
res9@tiYAxisString = "LST_Aqua_regrid(24,:,:)"
plot = gsn_csm_contour(wks1,LST_Aqua_regrid(24,:,:),res9)

; ================================== Tair ======================================
if diff .eq. "diff" then

  Tair     = new((/3,31,lat_sum,lon_sum/),float,default_val)
  Tair     = 0.

  filename1 = path(0)+"/cable_out_2019_SE_Aus.nc"
  f1 = addfile (filename1,"r")
  Tair_tmp   = f1->Tair(0:60*8,{A:B},{C:D})
  do event_num = 0,event_sum-1
    do day = HW_day_s(event_num),HW_day_e(event_num)
      ts = (day-1)*8
      ; to compare with Terra - 10am
      Tair(0,day-1,:,:) =  Tair_tmp(ts+3,:,:) + (Tair_tmp(ts+4,:,:) - Tair_tmp(ts+3,:,:))/3
      ; to compare with Aqua  - 2pm
      Tair(1,day-1,:,:) =  Tair_tmp(ts+4,:,:) + (Tair_tmp(ts+5,:,:) - Tair_tmp(ts+4,:,:))*(2/3)
      ; Max Tair
      Tair(2,day-1,:,:) =  dim_max_n_Wrap(Tair_tmp(ts:ts+7,:,:),0)
    end do
  end do
  delete(Tair_tmp)
  Tair = Tair - 273.15

  VegT(0,:,:,:,:) = VegT(0,:,:,:,:) - Tair(0:1,:,:,:)
  VegT(1,:,:,:,:) = VegT(1,:,:,:,:) - Tair(0:1,:,:,:)
  VegT(2,:,:,:,:) = VegT(2,:,:,:,:) - Tair(0:1,:,:,:)
  LST_Terra_regrid= LST_Terra_regrid - Tair(0,:,:,:)
  LST_Aqua_regrid = LST_Aqua_regrid  - Tair(1,:,:,:)
end if

; =============================== Diurnal Cycle ================================
; ---------------- Mask ----------------
; mask trees
filename1 = path(0)+"/cable_out_2019_SE_Aus.nc"
f1        = addfile (filename1,"r")

PFT_mask        = new((/3,lat_sum,lon_sum/),float,f1->iveg@_FillValue)
PFT_mask(0,:,:) = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, -1) ; PFT_mask@_FillValue)
PFT_mask(1,:,:) = where(  \
                        (f1->latitude({A:B},{C:D}) .ge. A_rgn_1)  .and. (f1->latitude({A:B},{C:D}) .le. B_rgn_1)   \
                  .and. (f1->longitude({A:B},{C:D}) .ge. C_rgn_1) .and. (f1->longitude({A:B},{C:D}) .le. D_rgn_1), \
                  PFT_mask(0,:,:),-1) ; PFT_mask@_FillValue)
PFT_mask(2,:,:) = where(  \
                        (f1->latitude({A:B},{C:D}) .ge. A_rgn_2)  .and. (f1->latitude({A:B},{C:D}) .le. B_rgn_2)   \
                  .and. (f1->longitude({A:B},{C:D}) .ge. C_rgn_2) .and. (f1->longitude({A:B},{C:D}) .le. D_rgn_2),  \
                  PFT_mask(0,:,:), -1); PFT_mask@_FillValue)

PFT_mask_multi = new((/9,3,lat_sum,lon_sum/),float,PFT_mask@_FillValue)
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
filename_hw = "/g/data/w35/mm3972/scripts/ehfheatwaves/HW_Event_Indicator_1970-2019.nc" ; summer HW only
fhw         = addfile(filename_hw,"r")

HW = new((/31,lat_sum,lon_sum/),float)
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

HW_mask_multi = new((/3,9,lat_sum,lon_sum/),float)
do z = 0, 8
  ; 15 Jan 2019
  HW_mask_multi(0,z,:,:) = HW(14,:,:)
  ; 16 Jan 2019
  HW_mask_multi(1,z,:,:) = HW(15,:,:)
  ; 25 Jan 2019
  HW_mask_multi(2,z,:,:) = HW(24,:,:)
end do
; delete(HW)

res9@tiYAxisString = "HW_mask_multi(0,0,:,:)"
plot = gsn_csm_contour(wks1,HW_mask_multi(0,0,:,:),res9)
res9@tiYAxisString = "HW_mask_multi(1,2,:,:)"
plot = gsn_csm_contour(wks1,HW_mask_multi(1,2,:,:),res9)
res9@tiYAxisString = "HW_mask_multi(2,4,:,:)"
plot = gsn_csm_contour(wks1,HW_mask_multi(2,4,:,:),res9)

; LST
LST_mask = new((/31,9,lat_sum,lon_sum/),float,-9999.)
LST_mask = 0.
do z = 0, 8
  LST_mask(:,z,:,:) = where(LST_Terra_regrid .ne. -9999. .and. LST_Aqua_regrid .ne. -9999., 1., 2.)
end do
res9@tiYAxisString = "LST_mask(14,0,:,:)"
plot = gsn_csm_contour(wks1,LST_mask(14,0,:,:),res9)
res9@tiYAxisString = "LST_mask(15,0,:,:)"
plot = gsn_csm_contour(wks1,LST_mask(15,0,:,:),res9)
res9@tiYAxisString = "LST_mask(24,0,:,:)"
plot = gsn_csm_contour(wks1,LST_mask(24,0,:,:),res9)

; ---------------- Diurnal Cycle ----------------
Fwsoil   = new((/case_sum,3,3,9/), float, default_val) ; case x day x region x time steps
deltaT   = new((/case_sum,3,3,9/), float, default_val)
Tair_ts  = new((/case_sum,3,3,9/), float, default_val)
Qair_ts  = new((/case_sum,3,3,9/), float, default_val)
VPD      = new((/case_sum,3,3,9/), float, default_val)
Fwsoil   = 0.
deltaT   = 0.
Tair_ts  = 0.
Qair_ts  = 0.
VPD      = 0.

do case_num = 0, case_sum-1
  print("case_num = " +case_num)
  filename1 = path(case_num)+"/cable_out_2019_SE_Aus.nc"
  f1 = addfile (filename1,"r")

  Tair_tmp   = f1->Tair(0:60*8,{A:B},{C:D})
  VegT_tmp   = f1->VegT(0:60*8,{A:B},{C:D})
  Qair_tmp   = f1->Qair(0:60*8,{A:B},{C:D})
  Fwsoil_tmp = f1->Fwsoil(0:60*8,{A:B},{C:D})

  do day = 0, dimsizes(HW_day_s) - 1
    ts = (HW_day_s(day)-1)*8
    do region = 0,2
      deltaT(case_num,day,region,:) = dim_avg_n_Wrap(where(\
                                              LST_mask(HW_day_s(day)-1,:,:,:) .eq. 1. \
                                              .and. PFT_mask_multi(:,region,:,:) .eq. 1 \
                                              .and. HW_mask_multi(day,:,:,:) .eq. 1, \
                                              VegT_tmp(ts:ts+8,:,:) - Tair_tmp(ts:ts+8,:,:), default_val) ,(/1,2/))
      Fwsoil(case_num,day,region,:) = dim_avg_n_Wrap(where(\
                                              LST_mask(HW_day_s(day)-1,:,:,:) .eq. 1. \
                                              .and. PFT_mask_multi(:,region,:,:) .eq. 1 \
                                              .and. HW_mask_multi(day,:,:,:) .eq. 1, \
                                              Fwsoil_tmp(ts:ts+8,:,:), default_val) ,(/1,2/))
      Tair_ts(case_num,day,region,:)= dim_avg_n_Wrap(where(\
                                              LST_mask(HW_day_s(day)-1,:,:,:) .eq. 1. \
                                              .and. PFT_mask_multi(:,region,:,:) .eq. 1 \
                                              .and. HW_mask_multi(day,:,:,:) .eq. 1, \
                                              Tair_tmp(ts:ts+8,:,:), default_val) ,(/1,2/))
      Qair_ts(case_num,day,region,:)= dim_avg_n_Wrap(where(\
                                              LST_mask(HW_day_s(day)-1,:,:,:) .eq. 1. \
                                              .and. PFT_mask_multi(:,region,:,:) .eq. 1 \
                                              .and. HW_mask_multi(day,:,:,:) .eq. 1, \
                                              Qair_tmp(ts:ts+8,:,:), default_val) ,(/1,2/))
    end do
  end do
  delete(filename1)
  delete(f1)
end do ; case_num

; ---------------- unit convert -----------------
Tair_ts  = Tair_ts -  273.15

; ---------------- Calculate VPD ----------------
Purf     = 1000.*100.

; saturation vapor pressure
Es = 100.0 * 6.112 * exp((17.67 * Tair_ts) / (243.5 + Tair_ts))

; vapor pressure
Ea = (Qair_ts * Purf) / (0.622 + (1.0 - 0.622) * Qair_ts)

VPD = (Es - Ea) * 0.001 ; in kPa
VPD = where(VPD .lt. 0.05, 0.05, VPD)

; ----------------- LST dots ------------------
LST_tree = new((/3,3,2/),float,-9999.) ; events x regions x 2 time steps
do region = 0,2
  LST_tree(0,region,0) = avg(where( PFT_mask(region,:,:) .eq. 1 .and. HW(14,:,:) .eq. 1\
                             .and. LST_Aqua_regrid(14,:,:) .ne. -9999., LST_Terra_regrid(14,:,:), -9999.) )
  LST_tree(1,region,0) = avg(where( PFT_mask(region,:,:) .eq. 1 .and. HW(15,:,:) .eq. 1\
                             .and. LST_Aqua_regrid(15,:,:) .ne. -9999., LST_Terra_regrid(15,:,:), -9999.) )
  LST_tree(2,region,0) = avg(where( PFT_mask(region,:,:) .eq. 1 .and. HW(24,:,:) .eq. 1\
                             .and. LST_Aqua_regrid(24,:,:) .ne. -9999., LST_Terra_regrid(24,:,:), -9999.) )
  LST_tree(0,region,1) = avg(where( PFT_mask(region,:,:) .eq. 1 .and. HW(14,:,:) .eq. 1\
                             .and. LST_Terra_regrid(14,:,:) .ne. -9999., LST_Aqua_regrid(14,:,:), -9999.) )
  LST_tree(1,region,1) = avg(where( PFT_mask(region,:,:) .eq. 1 .and. HW(15,:,:) .eq. 1\
                             .and. LST_Terra_regrid(15,:,:) .ne. -9999., LST_Aqua_regrid(15,:,:), -9999.) )
  LST_tree(2,region,1) = avg(where( PFT_mask(region,:,:) .eq. 1 .and. HW(24,:,:) .eq. 1\
                             .and. LST_Terra_regrid(24,:,:) .ne. -9999., LST_Aqua_regrid(24,:,:), -9999.) )
end do

; delete(res9)

; ================================= Plotting ===================================
pic = "./plots/"+pic_name
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
res@lbTitleFontHeightF     = 0.013
res@lbLabelFontHeightF     = 0.013

res@tmXBLabelFont           = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont           = 25
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF    = 0.015
res@tmYLLabelFontHeightF    = 0.015
; res@mpDataBaseVersion       = "HighRes"

res@lbLabelBarOn            = True
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

;**************** ͼ�������� *****************
pres                    = True
pres@wkWidth            = 2500
pres@wkHeight           = 2500

pres@gsnMaximize        = False
pres@gsnPanelLabelBar   = False
; pres@gsnPanelMainPosXF  = 0.1
; pres@gsnPanelMainPosYF  = 0.1
pres@cnLevelSelectionMode = "ExplicitLevels"
pres@lbBoxLinesOn       = True
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01
pres@gsnPanelYWhiteSpacePercent = 0.6
pres@gsnPanelXWhiteSpacePercent = 0.5


; --------- lines in diurnal cycle ---------
res1                    = True
res1@gsnDraw            = False
res1@gsnFrame           = False
res1@gsnMaximize        = False
res1@vpHeightF          = 0.32               ; change aspect ratio of plot
res1@vpWidthF 	        = 0.38
; res1@vpXF               = 0.08
; res1@vpYF               = 0.7
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
res1@pmLegendParallelPosF   = 0.25
res1@pmLegendOrthogonalPosF = -1.05
res1@pmLegendWidthF         = 0.10
res1@pmLegendHeightF        = 0.07

res1@tmXBMode            = "Explicit"
x                        = ispan(1,9,1)
res1@tmXBValues          = ispan(1,9,1)
res1@tmXBLabels          = tostring(ispan(0,24,3))
res1@trXMinF             = 1
res1@trXMaxF             = 9


; --------- dots in diurnal cycle ---------
res2 = res1
res2@xyMarkLineMode         = "Markers"
res2@xyMarkers              = (/6, 16/)
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
delete(res2@lgItemOrder)
res2@lgItemOrder            = (/ 1,0 /)   ; Reorder the legends

; -------------- gsn_panel ----------------
pres                    = True
pres@wkWidth            = 2500
pres@wkHeight           = 2500

pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False
pres@cnLevelSelectionMode = "ExplicitLevels"
pres@lbBoxLinesOn       = True
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01
pres@gsnPanelYWhiteSpacePercent = 0.7
pres@gsnPanelXWhiteSpacePercent = 0.5

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
plots = new(10,graphic)
plots2 = new(10,graphic)
; -------------- Tmax --------------
; 15 Jan
; res@cnLevels      = fspan( 20, 60, 41) ;deltaT_CTL_lvl
; res@cnFillPalette = "WhiteBlueGreenYellowRed"
; plots(0) = gsn_csm_contour_map(wks,Tair(2,14,:,:),res)
; ; 16 Jan
; plots(1) = gsn_csm_contour_map(wks,Tair(2,15,:,:),res)
; ; 25 Jan
; plots(2) = gsn_csm_contour_map(wks,Tair(2,24,:,:),res)
; delete(res@cnLevels)

; -------------- LST --------------
res@cnFillPalette = "precip3_16lev";"testcmap";"ViBlGrWhYeOrRe";"WhiteBlueGreenYellowRed"
res@cnLevels      = fspan(0,15, 16) ;deltaT_CTL_lvl

res@lbLabelBarOn = False ; "Always"
tmp = LST_Aqua_regrid(14,:,:)
tmp = where(PFT_mask(0,:,:) .eq. 1, tmp, tmp@_FillValue)
plots(0) = gsn_csm_contour_map(wks,tmp,res) ;LST_Aqua_regrid(14,:,:),res)
; plots(1) = gsn_csm_contour_map(wks,LST_Aqua_regrid(15,:,:),res)
delete(tmp)

res@lbLabelBarOn = True
tmp = LST_Aqua_regrid(24,:,:)
tmp = where(PFT_mask(0,:,:) .eq. 1, tmp, tmp@_FillValue)
plots(1) = gsn_csm_contour_map(wks,tmp,res) ;LST_Aqua_regrid(24,:,:),res)
; delete(res@cnLevels)
delete(tmp)

; -------------- deltaT --------------
; notes :
; VegT: case x 2 time steps x days x lat x lon

res@lbLabelBarOn = False ; "Always"
plots(2) = gsn_csm_contour_map(wks,VegT(0,1,14,:,:),res)
; plots(4) = gsn_csm_contour_map(wks,VegT(0,1,15,:,:),res)
res@lbLabelBarOn = True ; "Always"
plots(3) = gsn_csm_contour_map(wks,VegT(0,1,24,:,:),res)
delete(res@cnLevels)

; -------------- deltaT GW - FD --------------
res@cnLevels      = (/-5,-4.5,-4,-3.5,-3.,-2.5,-2.,-1.5,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5./)
                   ;fspan(-4,4, 17)
res@cnFillPalette = "BlAqGrWh2YeOrReVi22";"BlAqGrWh2YeOrReVi22";"testcmap";"ViBlGrWhYeOrRe";"WhiteBlueGreenYellowRed"
tmp               = VegT(0,1,:,:,:)
tmp               = tmp - VegT(1,1,:,:,:)
; temp = VegT(0,1,0,:,:)
; temp = where(PFT_mask(1,:,:) .eq. 1 .and. HW(14,:,:) .eq. 1 .and. LST_mask(14,0,:,:) .eq. 1,\
;       tmp(14,:,:), tmp@_FillValue)
res@lbLabelBarOn = False ; "Always"
plots(4) = gsn_csm_contour_map(wks,tmp(14,:,:),res)
; delete(temp)

; plots(7) = gsn_csm_contour_map(wks,VegT(0,1,15,:,:),res)
res@lbLabelBarOn = True ; "Always"
; temp = where(PFT_mask(1,:,:) .eq. 1 .and. HW(24,:,:) .eq. 1 .and. LST_mask(24,0,:,:) .eq. 1,\
;       tmp(24,:,:), tmp@_FillValue)
plots(5) = gsn_csm_contour_map(wks,tmp(24,:,:),res)
; delete(res@cnLevels)
delete(tmp)
; delete(temp)

; -------------- deltaT DR - GW --------------
res@cnFillPalette = "BlAqGrWh2YeOrReVi22"
tmp               = VegT(2,1,:,:,:)
tmp               = tmp - VegT(0,1,:,:,:)
res@lbLabelBarOn = False ; "Always"
plots2(0) = gsn_csm_contour_map(wks,tmp(14,:,:),res)
; plots(7) = gsn_csm_contour_map(wks,VegT(0,1,15,:,:),res)
res@lbLabelBarOn = True ; "Always"
plots2(1) = gsn_csm_contour_map(wks,tmp(24,:,:),res)
delete(res@cnLevels)
delete(tmp)

; -------------- diurnal cycle --------------
; notes:
; deltaT : case x day (0,1,2) x region x time steps
; LST_tree : events x regions x 2 time steps
offset_num = 6;8
; plots = new(4,graphic)
res1@pmLegendDisplayMode = "Always"
res1@trYMinF             = -2.2
res1@trYMaxF             = 11.
y                        = (/4.33, 5.67/) ; 10 am & 2 pm
res1@tiYAxisString       = "~F5~D~F22~T (~S~o~N~C)"

res1@tiXAxisString       = "hr"
; region 1
print("region1 15th")
print(deltaT(:,0,1,:))
print(LST_tree(0,1,:))
plots(0+offset_num) = gsn_csm_xy(wks,x,deltaT(:,0,1,:),res1)
plot1    = gsn_csm_xy(wks,y,LST_tree(0,1,:),res2)
overlay(plots(0+offset_num),plot1)
delete(res1@tiYAxisString)

res1@pmLegendDisplayMode = "Never"
; plots(10) = gsn_csm_xy(wks,x,deltaT(:,1,:),res1)
; plot2     = gsn_csm_xy(wks,y,LST_tree(1,:),res2)
; overlay(plots(10),plot2)
;
; print("region1 25th")
; print(deltaT(:,2,1,:))
; print(LST_tree(2,1,:))
; plots(1+offset_num) = gsn_csm_xy(wks,x,deltaT(:,2,1,:),res1)
; plot2     = gsn_csm_xy(wks,y,LST_tree(2,1,:),res2)
; overlay(plots(1+offset_num),plot2)
;
; ; region 2
; res1@tiYAxisString       = "~F5~D~F22~T (~S~o~N~C)"
; res1@tiXAxisString       = "hr"
; print("region2 15th")
; print(deltaT(:,0,2,:))
; print(LST_tree(0,2,:))
; plots(2+offset_num) = gsn_csm_xy(wks,x,deltaT(:,0,2,:),res1)
; plot3    = gsn_csm_xy(wks,y,LST_tree(0,2,:),res2)
; overlay(plots(2+offset_num),plot3)
; delete(res1@tiYAxisString)

; res1@pmLegendDisplayMode = "Never"
; plots(10) = gsn_csm_xy(wks,x,deltaT(:,1,:),res1)
; plot2     = gsn_csm_xy(wks,y,LST_tree(1,:),res2)
; overlay(plots(10),plot2)
res1@tiXAxisString       = "hr"
print("region2 15th")
print(deltaT(:,2,2,:))
print(LST_tree(2,2,:))
plots(1+offset_num) = gsn_csm_xy(wks,x,deltaT(:,2,2,:),res1) ; 3+offset_num
plot4     = gsn_csm_xy(wks,y,LST_tree(2,2,:),res2)
overlay(plots(1+offset_num),plot4)

; -------------- VPD --------------
; res1@tiYAxisString       = "VPD (kPa)"
; res1@trYMinF             = 0.
; res1@trYMaxF             = 7.
; plots(12) = gsn_csm_xy(wks,x,VPD(:,0,:),res1)
; delete(res1@tiYAxisString)
; plots(13) = gsn_csm_xy(wks,x,VPD(:,1,:),res1)
; plots(14) = gsn_csm_xy(wks,x,VPD(:,2,:),res1)
;
; res1@tiYAxisString       = "~F5~b"
; res1@trYMinF             = 0.
; res1@trYMaxF             = 1.
; plots(15) = gsn_csm_xy(wks,x,Fwsoil(:,0,:),res1)
; delete(res1@tiYAxisString)
; plots(16) = gsn_csm_xy(wks,x,Fwsoil(:,1,:),res1)
; plots(17) = gsn_csm_xy(wks,x,Fwsoil(:,2,:),res1)

; ________________ add shading ________________
poly             = new(6,graphic)
gres             = True
gres@gsEdgesOn   = True
gres@gsEdgeColor = "red"
gres@gsLineColor = "red"
gres@gsLineThicknessF = 3.0
; gres@gsFillColor = "cornsilk" ;"rosybrown" ;"azure"
; gres@gsFillOpacityF = 0.0

xmin = (/C_rgn_1,C_rgn_2/)
xmax = (/D_rgn_1,D_rgn_2/)
ymin = (/A_rgn_1,A_rgn_2/)
ymax = (/B_rgn_1,B_rgn_2/)

xbars = (/ (/xmin(0), xmax(0), xmax(0), xmin(0), xmin(0)/), \
            (/xmin(1), xmax(1), xmax(1), xmin(1), xmin(1)/) /)

ybars = (/ (/ymin(0),ymin(0),ymax(0),ymax(0),ymin(0)/), \
           (/ymin(1),ymin(1),ymax(1),ymax(1),ymin(1)/)  /)

do counter = 0, 5
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

; _____________ DR-GW _______________
poly2 = new(2,graphic)
do cnt = 0, 1
  txid1   = gsn_create_text(wks, order_name(cnt), txres)
  annoid1 = gsn_add_annotation(plots2(cnt), txid1, amres)
end do

do counter = 0, 1
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

; ------------------ Check LST missing value -----------------
; plots = new(2,graphic)
; LST_Terra_regrid = where(HW .eq. 1, LST_Terra_regrid,LST_Terra_regrid@_FillValue)
; LST_Aqua_regrid  = where(HW .eq. 1, LST_Aqua_regrid,LST_Aqua_regrid@_FillValue)
; do day = 0, 30
;   LST_Terra_regrid(day,:,:) = where(PFT_mask .eq. 1, LST_Terra_regrid(day,:,:),LST_Terra_regrid@_FillValue)
;   LST_Aqua_regrid(day,:,:)  = where(PFT_mask .eq. 1, LST_Aqua_regrid(day,:,:),LST_Aqua_regrid@_FillValue)
; end do
; do day = 0,30
;   print("VegT")
;   res@tiXAxisString = "day = " +(day+1)
;   plots(0) = gsn_csm_contour_map(wks,LST_Terra_regrid(day,:,:),res)
;   plots(1) = gsn_csm_contour_map(wks,LST_Aqua_regrid(day,:,:),res)
;   gsn_panel(wks,(/plots/),(/1,2/),pres)
; end do

; ------------------ Check Jan temp ------------------
; do day = 0,30
;   print("VegT")
;   plots(0) = gsn_csm_contour_map(wks,VegT(0,0,day,:,:),res)
;   plots(1) = gsn_csm_contour_map(wks,VegT(1,0,day,:,:),res)
;   plots(2) = gsn_csm_contour_map(wks,VegT(2,0,day,:,:),res)
;   print("LST")
;   plots(3) = gsn_csm_contour_map(wks,LST_Terra_regrid(day,:,:),res)
;
;   print("VegT")
;   plots(4) = gsn_csm_contour_map(wks,VegT(0,1,day,:,:),res)
;   plots(5) = gsn_csm_contour_map(wks,VegT(1,1,day,:,:),res)
;   plots(6) = gsn_csm_contour_map(wks,VegT(2,1,day,:,:),res)
;   print("LST")
;   plots(7) = gsn_csm_contour_map(wks,LST_Aqua_regrid(day,:,:),res)
;
;   do cnt = 0, 7
;     txid1   = gsn_create_text(wks, order_name(cnt), txres)
;     annoid1 = gsn_add_annotation(plots(cnt), txid1, amres)
;   end do
;   gsn_panel(wks,(/plots/),(/2,4/),pres)
; end do


end
