;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

; ================================== Options ===================================
year_s    = 2017
year_e    = 2019
Soil_dryness = True
PFT       = "all" ; "tree" ;"all"
plot_type = "all"   ;"water_bal" ;"all"
scale     = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
heatwave  = "hw-90th" ; "hw-90th" ; "all-day"
case_name = (/"CTL","FREE_DRAIN"/);,"PUMPx10","PUMPx100"/)

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

pic_name  = "HW_ET_deltaT_"+heatwave+"_PFT-"+PFT+"_"+scale+"_spinup200yr_"+year_s+"-"+year_e

; =================================== Path =====================================
path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"+case_name(case_num)+"/outputs-spinup200yr"
end do

; =========================== PFT & Spitial Scale ==============================
; set PFT
if PFT .eq. "all" then
  iveg_s = 0
  iveg_e = 20
else if PFT .eq. "tree" then
  iveg_s = 1
  iveg_e = 4
end if
end if

if scale .eq. "Global" then
  A = -90.
  B = 90.
  C = 0.
  D = 360.
else if scale .eq. "SE-AU" then
  A = -40.
  B = -28.
  C = 140.
  D = 154.
else if scale .eq. "AU" then
  A = -44.
  B = -10.
  C = 112.
  D = 154.
end if
end if
end if

; ================================== Start =====================================
filename = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f        = addfile (filename,"r")
lat_sum  = dimsizes(f->Rainf(0,{A:B},0))
lon_sum  = dimsizes(f->Rainf(0,0,{C:D}))
def_value = f->Rainf@_FillValue

; PFT mask
PFT_mask  = where((f->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
PFT_mask_multi = new((/day_sum,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
do z = 0, day_sum-1
  PFT_mask_multi(z,:,:) = PFT_mask
end do
delete(PFT_mask)
delete(f)

Rnet     = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
TVeg     = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
Fwsoil   = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
Tair     = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
VegT     = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
SMsurf   = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
SM64cm   = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
; SoilMoist= new((/case_sum,day_sum,6,lat_sum,lon_sum/),float,def_value)
WatTable = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)

Rnet     = 0.
TVeg     = 0.
Fwsoil   = 0.
Tair     = 0.
VegT     = 0.
SMsurf   = 0.
SM64cm   = 0.
; SoilMoist= 0.
WatTable = 0.

; _______ READ CABLE ________
day_step_s = 0
do year = year_s,year_e
  print("year = " +year)

    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    yr_tot_day = sum(dom)
    day_step_e = day_step_s + yr_tot_day - 1
    do case_num = 0,case_sum-1
      print("case_num = " +case_num)
      filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
      f1 = addfile (filename1,"r")
      TVeg(case_num,day_step_s:day_step_e,:,:)     = f1->TVeg(:,{A:B},{C:D})*3600.*24.
      Fwsoil(case_num,day_step_s:day_step_e,:,:)   = f1->Fwsoil(:,{A:B},{C:D})
      WatTable(case_num,day_step_s:day_step_e,:,:) = f1->WatTable(:,{A:B},{C:D})
      Tair(case_num,day_step_s:day_step_e,:,:)     = f1->Tair(:,{A:B},{C:D})-273.15
      VegT(case_num,day_step_s:day_step_e,:,:)     = f1->VegT(:,{A:B},{C:D})-273.15
      SMsurf(case_num,day_step_s:day_step_e,:,:)   = ( f1->SoilMoist(:,0,{A:B},{C:D})*0.022 \
                                                     + f1->SoilMoist(:,1,{A:B},{C:D})*0.058 \
                                                     + f1->SoilMoist(:,2,{A:B},{C:D})*(0.1 - 0.022 - 0.058)) \
                                                     / 0.1
      SM64cm(case_num,day_step_s:day_step_e,:,:)   = ( f1->SoilMoist(:,0,{A:B},{C:D})*0.022  \
                                                      + f1->SoilMoist(:,1,{A:B},{C:D})*0.058 \
                                                      + f1->SoilMoist(:,2,{A:B},{C:D})*0.154 \
                                                      + f1->SoilMoist(:,3,{A:B},{C:D})*0.409)\
                                                      / 0.643
      Rnet(case_num,day_step_s:day_step_e,:,:)     = f1->Rnet(:,{A:B},{C:D})
      ; SoilMoist(case_num,day_step_s:day_step_e,:,:,:)= f1->SoilMoist(:,:,{A:B},{C:D})
      delete(f1)
    end do ; case_num
    if heatwave .ne. "all-day" then
       day_step_s = day_step_e + 1
    end if
end do ; year

; mask PFT
print("mask PFT")
do case_num = 0,case_sum-1
  TVeg(case_num,:,:,:)     = where(PFT_mask_multi .eq. 1, TVeg(case_num,:,:,:), def_value)
  Fwsoil(case_num,:,:,:)   = where(PFT_mask_multi .eq. 1, Fwsoil(case_num,:,:,:), def_value)
  WatTable(case_num,:,:,:) = where(PFT_mask_multi .eq. 1, WatTable(case_num,:,:,:), def_value)
  Tair(case_num,:,:,:)     = where(PFT_mask_multi .eq. 1, Tair(case_num,:,:,:), def_value)
  VegT(case_num,:,:,:)     = where(PFT_mask_multi .eq. 1, VegT(case_num,:,:,:), def_value)
  SMsurf(case_num,:,:,:)   = where(PFT_mask_multi .eq. 1, SMsurf(case_num,:,:,:), def_value)
  SM64cm(case_num,:,:,:)   = where(PFT_mask_multi .eq. 1, SM64cm(case_num,:,:,:), def_value)
  Rnet(case_num,:,:,:)     = where(PFT_mask_multi .eq. 1, Rnet(case_num,:,:,:), def_value)
  ; do soil = 0, 5
    ; SoilMoist(case_num,:,soil,:,:)= where(PFT_mask_multi .eq. 1, SoilMoist(case_num,:,soil,:,:), def_value)
  ; end do
end do
delete(PFT_mask_multi)

; mask HW
print("mask HW")
if heatwave .eq. "hw-90th" then
  filename_hw = "/g/data/w35/mm3972/scripts/ehfheatwaves/HW_Event_Indicator_1970-2019.nc" ; summer HW only
  fhw         = addfile(filename_hw,"r")

  HW = new((/day_sum,lat_sum,lon_sum/),float)
  day_start = 0
  do year = 1970, year_s - 1
    if ( mod(year,4) .eq. 0) then
       day_start = day_start + 366
    else
       day_start = day_start + 365
    end if
  end do
  HW = fhw->event(day_start:day_start+day_sum-1,:,:)*1.
  delete(filename_hw)
  delete(fhw)

  do case_num = 0,case_sum-1
    print(case_name(case_num))
    TVeg(case_num,:,:,:)      = where(HW .eq. 1, TVeg(case_num,:,:,:), def_value)
    Fwsoil(case_num,:,:,:)    = where(HW .eq. 1, Fwsoil(case_num,:,:,:), def_value)
    WatTable(case_num,:,:,:)  = where(HW .eq. 1, WatTable(case_num,:,:,:), def_value)
    Tair(case_num,:,:,:)      = where(HW .eq. 1, Tair(case_num,:,:,:), def_value)
    VegT(case_num,:,:,:)      = where(HW .eq. 1, VegT(case_num,:,:,:), def_value)
    SMsurf(case_num,:,:,:)    = where(HW .eq. 1, SMsurf(case_num,:,:,:), def_value)
    SM64cm(case_num,:,:,:)    = where(HW .eq. 1, SM64cm(case_num,:,:,:), def_value)
    Rnet(case_num,:,:,:)      = where(HW .eq. 1, Rnet(case_num,:,:,:), def_value)
  end do
  delete(HW)
end if

TVeg_1d      = new((/case_sum,day_sum/), float, def_value)
Fwsoil_1d    = new((/case_sum,day_sum/), float, def_value)
WatTable_1d  = new((/case_sum,day_sum/), float, def_value)
Tair_1d      = new((/case_sum,day_sum/), float, def_value)
VegT_1d      = new((/case_sum,day_sum/), float, def_value)
SMsurf_1d    = new((/case_sum,day_sum/), float, def_value)
SM64cm_1d    = new((/case_sum,day_sum/), float, def_value)
Rnet_1d      = new((/case_sum,day_sum/), float, def_value)

TVeg_1d      = dim_avg_n_Wrap(TVeg,(/2,3/))
Fwsoil_1d    = dim_avg_n_Wrap(Fwsoil,(/2,3/))
WatTable_1d  = dim_avg_n_Wrap(WatTable,(/2,3/))
Tair_1d      = dim_avg_n_Wrap(Tair,(/2,3/))
VegT_1d      = dim_avg_n_Wrap(VegT,(/2,3/))
SMsurf_1d    = dim_avg_n_Wrap(SMsurf,(/2,3/))
SM64cm_1d    = dim_avg_n_Wrap(SM64cm,(/2,3/))
Rnet_1d      = dim_avg_n_Wrap(Rnet,(/2,3/))

delete(WatTable)
delete(SMsurf)
delete(SM64cm)
delete(Rnet)

if Soil_dryness then
  print("Soil_dryness")
  TVeg_1d_SM      = new((/3,case_sum,day_sum/), float, def_value)
  Tair_1d_SM      = new((/3,case_sum,day_sum/), float, def_value)
  VegT_1d_SM      = new((/3,case_sum,day_sum/), float, def_value)
  TVeg_1d_SM(0,:,:) = where(Fwsoil_1d .ge. 0.8, TVeg_1d, def_value)
  TVeg_1d_SM(1,:,:) = where(Fwsoil_1d .gt. 0.3 .and. Fwsoil_1d .lt. 0.8, TVeg_1d, def_value)
  TVeg_1d_SM(2,:,:) = where(Fwsoil_1d .le. 0.3, TVeg_1d, def_value)
  delete(TVeg)
  Tair_1d_SM(0,:,:) = where(Fwsoil_1d .ge. 0.8, Tair_1d, def_value)
  Tair_1d_SM(1,:,:) = where(Fwsoil_1d .gt. 0.3 .and. Fwsoil_1d .lt. 0.8, Tair_1d, def_value)
  Tair_1d_SM(2,:,:) = where(Fwsoil_1d .le. 0.3, Tair_1d, def_value)
  delete(Tair)
  VegT_1d_SM(0,:,:) = where(Fwsoil_1d .ge. 0.8, VegT_1d, def_value)
  VegT_1d_SM(1,:,:) = where(Fwsoil_1d .gt. 0.3 .and. Fwsoil_1d .lt. 0.8, VegT_1d, def_value)
  VegT_1d_SM(2,:,:) = where(Fwsoil_1d .le. 0.3, VegT_1d, def_value)
  delete(VegT)
end if
delete(Fwsoil)

; length       = dimsizes(ndtooned(TVeg(0,:,:,:)))
; print(length)
; TVeg_1d      = new((/case_sum,length/), float, def_value)
; Fwsoil_1d    = new((/case_sum,length/), float, def_value)
; WatTable_1d  = new((/case_sum,length/), float, def_value)
; Tair_1d      = new((/case_sum,length/), float, def_value)
; VegT_1d      = new((/case_sum,length/), float, def_value)
; SMsurf_1d    = new((/case_sum,length/), float, def_value)
; SM64cm_1d    = new((/case_sum,length/), float, def_value)
; Rnet_1d      = new((/case_sum,length/), float, def_value)
; ; SoilMoist_1d = new((/case_sum,6,length/), float, def_value)
;
; do case_num = 0,case_sum-1
;   print(case_name(case_num))
;   TVeg_1d(case_num,:)      = dimsizes(ndtooned(TVeg(case_num,:,:,:)))
;   Fwsoil_1d(case_num,:)    = dimsizes(ndtooned(Fwsoil(case_num,:,:,:)))
;   WatTable_1d(case_num,:)  = dimsizes(ndtooned(WatTable(case_num,:,:,:)))
;   Tair_1d(case_num,:)      = dimsizes(ndtooned(Tair(case_num,:,:,:)))
;   VegT_1d(case_num,:)      = dimsizes(ndtooned(VegT(case_num,:,:,:)))
;   SMsurf_1d(case_num,:)    = dimsizes(ndtooned(SMsurf(case_num,:,:,:)))
;   SM64cm_1d(case_num,:)    = dimsizes(ndtooned(SM64cm(case_num,:,:,:)))
;   Rnet_1d(case_num,:)      = dimsizes(ndtooned(Rnet(case_num,:,:,:)))
;   ; do soil = 0,5
;     ; SoilMoist_1d(case_num,soil,:) = dimsizes(ndtooned(SoilMoist(case_num,:,soil,:,:)))
;   ; end do
; end do
;
; delete(WatTable)
; delete(SMsurf)
; delete(SM64cm)
; delete(Rnet)
; ; delete(SoilMoist)

; if Soil_dryness then
;   print("Soil_dryness")
;   TVeg_1d_SM      = new((/3,case_sum,length/), float, def_value)
;   Tair_1d_SM      = new((/3,case_sum,length/), float, def_value)
;   VegT_1d_SM      = new((/3,case_sum,length/), float, def_value)
;   TVeg_1d_SM(0,:,:) = where(Fwsoil_1d .ge. 0.8, TVeg_1d, def_value)
;   TVeg_1d_SM(1,:,:) = where(Fwsoil_1d .gt. 0.3 .and. Fwsoil_1d .lt. 0.8, TVeg_1d, def_value)
;   TVeg_1d_SM(2,:,:) = where(Fwsoil_1d .le. 0.3, TVeg_1d, def_value)
;   delete(TVeg)
;   Tair_1d_SM(0,:,:) = where(Fwsoil_1d .ge. 0.8, Tair_1d, def_value)
;   Tair_1d_SM(1,:,:) = where(Fwsoil_1d .gt. 0.3 .and. Fwsoil_1d .lt. 0.8, Tair_1d, def_value)
;   Tair_1d_SM(2,:,:) = where(Fwsoil_1d .le. 0.3, Tair_1d, def_value)
;   delete(Tair)
;   VegT_1d_SM(0,:,:) = where(Fwsoil_1d .ge. 0.8, VegT_1d, def_value)
;   VegT_1d_SM(1,:,:) = where(Fwsoil_1d .gt. 0.3 .and. Fwsoil_1d .lt. 0.8, VegT_1d, def_value)
;   VegT_1d_SM(2,:,:) = where(Fwsoil_1d .le. 0.3, VegT_1d, def_value)
;   delete(VegT)
; end if
; delete(Fwsoil)
; ================================= Plotting ===================================
wks                      = gsn_open_wks("pdf","./plots/"+pic_name)
res                      = True                   ; plot mods desired
res@gsnDraw              = False  ; Don't draw plots
res@gsnFrame             = False  ; ����ҳ
res@gsnMaximize          = True                   ; maximize plot in frame
res@vpHeightF            = 0.90
res@vpWidthF             = 0.90
res@vpXF                 = 0.08
res@vpYF                 = 0.98
res@pmLegendWidthF       = 0.2                     ; resize legend width
res@pmLegendHeightF      = 0.2                     ; and height
res@tiXAxisFontHeightF   = 0.05
res@tiYAxisFontHeightF   = 0.05
res@tmXBLabelFontHeightF = 0.05
res@tmYLLabelFontHeightF = 0.05
res@gsnPaperOrientation = "portrait";"portrait";"landscape" ;"auto"

res@tmXBMode             = "Explicit"
res@tmYLMode             = "Explicit"
res@xyMarkLineModes      = (/"Markers","Markers","Markers","Markers"/) ;,"Lines","Lines"/)  ; choose which have markers
res@xyMarkers            = 16                     ; choose type of marker

res@xyMarkerColors       = (/ "red","blue","green","orange" /)
res@xyMarkerSizeF        = (/ 0.005 ,0.005, 0.005 ,0.005 /)

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

plots            = new(6,graphic)

; TVeg vs SM64cm
print("TVeg vs SM64cm")
res@tiXAxisString       = "SM64cm (m3/m3)"
res@tiYAxisString       = "TVeg (mm/d)"
res@trXMinF             = 0.
res@trXMaxF             = 0.4
res@trYMinF             = 0.
res@trYMaxF             = 5.
plots(0)                = gsn_csm_xy (wks,SM64cm_1d,TVeg_1d,res)        ; create plot

; TVeg vs SMsurf
print("TVeg vs SMsurf")
res@tiXAxisString       = "SMsurf (m3/m3)"
res@tiYAxisString       = "TVeg (mm/d)"
res@trXMinF             = 0.
res@trXMaxF             = 0.4
res@trYMinF             = 0.
res@trYMaxF             = 5.
plots(1)                = gsn_csm_xy (wks,SMsurf_1d,TVeg_1d,res)        ; create plot

; TVeg vs Tair
print("TVeg vs Tair")
res@tiXAxisString       = "Tair (m3/m3)"
res@tiYAxisString       = "TVeg (mm/d)"
res@trXMinF             = 10.
res@trXMaxF             = 35.
res@trYMinF             = 0.
res@trYMaxF             = 5.
plots(2)                = gsn_csm_xy (wks,Tair_1d,TVeg_1d,res)        ; create plot

; TVeg vs Rnet
print("TVeg vs Rnet")
res@tiXAxisString       = "Rnet (W/m2)"
res@tiYAxisString       = "TVeg (mm/d)"
res@trXMinF             = -50.
res@trXMaxF             = 300.
res@trYMinF             = 0.
res@trYMaxF             = 5.
plots(3)                = gsn_csm_xy (wks,Rnet_1d,TVeg_1d,res)        ; create plot

; TVeg vs Fwsoil
print("TVeg vs Fwsoil")
res@tiXAxisString       = "Fwsoil (-)"
res@tiYAxisString       = "TVeg (mm/d)"
res@trXMinF             = 0.
res@trXMaxF             = 1.
res@trYMinF             = 0.
res@trYMaxF             = 5.
plots(4)                = gsn_csm_xy (wks,Fwsoil_1d,TVeg_1d,res)        ; create plot

; TVeg vs WatTable
print("TVeg vs WatTable")
res@tiXAxisString       = "WatTable (m)"
res@tiYAxisString       = "TVeg (mm/d)"
res@trXMinF             = 0.
res@trXMaxF             = 20.
res@trYMinF             = 0.
res@trYMaxF             = 5.
plots(5)                = gsn_csm_xy (wks,WatTable_1d,TVeg_1d,res)        ; create plot

gsn_panel(wks,(/plots/),(/2,3/),pres)
delete(plots)

plots       = new(3,graphic)
; deltaT vs TVeg
print("deltaT vs TVeg")
res@tiXAxisString       = "TVeg (mm/d)"
res@tiYAxisString       = "VegT - Tair (C)"
res@trXMinF             = 0.
res@trXMaxF             = 5.
res@trYMinF             = -2.
res@trYMaxF             = 6.
plots(0)                = gsn_csm_xy (wks,TVeg_1d,VegT_1d - Tair_1d,res)        ; create plot

; deltaT vs Fwsoil
print("deltaT vs Fwsoil")
res@tiXAxisString       = "Fwsoil (-)"
res@tiYAxisString       = "VegT - Tair (C)"
res@trXMinF             = 0.
res@trXMaxF             = 1.
res@trYMinF             = -2.
res@trYMaxF             = 6.
plots(1)                = gsn_csm_xy (wks,Fwsoil_1d,VegT_1d - Tair_1d,res)        ; create plot

; deltaT vs SM64cm
print("deltaT vs SM64cm")
res@tiXAxisString       = "SM64cm (m3/m3)"
res@tiYAxisString       = "VegT - Tair (C)"
res@trXMinF             = 0.
res@trXMaxF             = 0.4
res@trYMinF             = -2.
res@trYMaxF             = 6.
plots(2)                 = gsn_csm_xy (wks,SM64cm_1d,VegT_1d - Tair_1d,res)        ; create plot

gsn_panel(wks,(/plots/),(/1,3/),pres)
delete(plots)


if Soil_dryness then
  plots            = new(6,graphic)
  ; VegT vs Tair
  res@tiMainString        = "Wet Soil (Fwsoil > 0.7)"
  res@tiXAxisString       = "Tair (C)"
  res@tiYAxisString       = "VegT (C)"
  res@trXMinF             = 10.
  res@trXMaxF             = 40.
  res@trYMinF             = 10.
  res@trYMaxF             = 40.
  plots(0)                 = gsn_csm_xy (wks,Tair_1d_SM(0,:,:),VegT_1d_SM(0,:,:),res)        ; create plot

  ; VegT vs Tair
  res@tiMainString        = "Wet Soil ( 0.3 < Fwsoil < 0.7)"
  res@tiXAxisString       = "Tair (C)"
  res@tiYAxisString       = "VegT (C)"
  res@trXMinF             = 10.
  res@trXMaxF             = 40.
  res@trYMinF             = 10.
  res@trYMaxF             = 40.
  plots(1)                 = gsn_csm_xy (wks,Tair_1d_SM(1,:,:),VegT_1d_SM(1,:,:),res)        ; create plot

  ; VegT vs Tair
  res@tiMainString        = "Wet Soil (Fwsoil < 0.3)"
  res@tiXAxisString       = "Tair (C)"
  res@tiYAxisString       = "VegT (C)"
  res@trXMinF             = 10.
  res@trXMaxF             = 40.
  res@trYMinF             = 10.
  res@trYMaxF             = 40.
  plots(2)                 = gsn_csm_xy (wks,Tair_1d_SM(2,:,:),VegT_1d_SM(2,:,:),res)        ; create plot

  ; VegT vs TVeg
  res@tiMainString        = "Wet Soil (Fwsoil > 0.7)"
  res@tiXAxisString       = "TVeg (mm/d)"
  res@tiYAxisString       = "VegT (C)"
  res@trXMinF             = 0.
  res@trXMaxF             = 5.
  res@trYMinF             = 10.
  res@trYMaxF             = 40.
  plots(3)                 = gsn_csm_xy (wks,TVeg_1d_SM(0,:,:),VegT_1d_SM(0,:,:),res)        ; create plot

  ; VegT vs TVeg
  res@tiMainString        = "Wet Soil ( 0.3 < Fwsoil < 0.7)"
  res@tiXAxisString       = "TVeg (mm/d)"
  res@tiYAxisString       = "VegT (C)"
  res@trXMinF             = 0.
  res@trXMaxF             = 5.
  res@trYMinF             = 10.
  res@trYMaxF             = 40.
  plots(4)                 = gsn_csm_xy (wks,TVeg_1d_SM(1,:,:),VegT_1d_SM(1,:,:),res)        ; create plot

  ; VegT vs TVeg
  res@tiMainString        = "Wet Soil (Fwsoil < 0.3)"
  res@tiXAxisString       = "TVeg (mm/d)"
  res@tiYAxisString       = "VegT (C)"
  res@trXMinF             = 0.
  res@trXMaxF             = 5.
  res@trYMinF             = 10.
  res@trYMaxF             = 40.
  plots(5)                 = gsn_csm_xy (wks,TVeg_1d_SM(2,:,:),VegT_1d_SM(2,:,:),res)        ; create plot

  gsn_panel(wks,(/plots/),(/2,3/),pres)
  delete(plots)
end if
end
