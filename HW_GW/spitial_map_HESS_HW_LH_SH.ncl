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

; ================================= Predefine ==================================
year_s   = 2000
year_e   = 2019
row_sum  = 4
col_sum  = 4
heatwave = "hw-90th" ; "hw-90th" ; "all-day"
scale    = "SE-AU"   ; "45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
case_name= (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr", \
             "CTL-30x3+30yr-DeepRoot/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr-DeepRoot"/)
path_mum = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
message  = "CTL_FD_DR_FD+DR-FD"
pic      = "./plots/Spitial_HESS_LH_SH_"+heatwave+"_"+message

year_sum = year_e-year_s+1
day_sum  = 0

do year = year_s,year_e
  if ( mod(year,4) .eq. 0) then
     day_sum = day_sum + 366
  else
     day_sum = day_sum + 365
  end if
end do

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

; =================================== Path =====================================
case_sum    = dimsizes(case_name)
path        = new(case_sum,"string")
do case_num = 0, case_sum -1
   path(case_num) = path_mum+case_name(case_num);+"/outputs"
end do

path_LST= "/g/data/w35/mm3972/data/MODIS_LST";

order_name  = (/"(a)","(b)","(c)","(d)",\
                "(e)","(f)","(g)","(h)",\
                "(i)","(j)","(k)","(l)",\
                "(m)","(n)","(o)","(p)"/)

;            all years,  2001-2009,  2017-2019,   2010-2011
;             2000-1-1,   2001-1-1,   2017-1-1,    2010-1-1
ts_s    = (/         0,        366,       6210,       3653/)
;           2019-12-31, 2009-12-31, 2019-12-31,  2011-12-31
ts_e    = (/ day_sum-1,       3652,       7304,       4382/)


; =================================== START ====================================
filename  = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f         = addfile(filename,"r")

lat_sum   = dimsizes(f->Rainf(0,{A:B},0))
lon_sum   = dimsizes(f->Rainf(0,0,{C:D}))
lat_size  = dimsizes(f->y)
lon_size  = dimsizes(f->x)
CABLE_lat = f->y({A:B})
CABLE_lon = f->x({C:D})
default_val = -1.e+33
delete(f)

Qle       = new((/day_sum, lat_sum, lon_sum/),float, default_val)
Qh        = new((/day_sum, lat_sum, lon_sum/),float, default_val)
Fwsoil    = new((/day_sum, lat_sum, lon_sum/),float, default_val)
deltaT    = new((/day_sum, lat_sum, lon_sum/),float, default_val)

Qle       = 0.
Qh        = 0.
Fwsoil    = 0.
deltaT    = 0.

time_period = 4
Qle_avg    = new((/case_sum, time_period, lat_sum, lon_sum/), float, default_val)
Qh_avg     = new((/case_sum, time_period, lat_sum, lon_sum/), float, default_val)
Fwsoil_avg = new((/case_sum, time_period, lat_sum, lon_sum/), float, default_val)
deltaT_avg = new((/case_sum, time_period, lat_sum, lon_sum/), float, default_val)

Qle_avg!0  = "case"
Qle_avg!1  = "time_period"
Qle_avg!2  = "lat"
Qle_avg!3  = "lon"
Qle_avg&lat= CABLE_lat
Qle_avg&lon= CABLE_lon

Qh_avg!0  = "case"
Qh_avg!1  = "time_period"
Qh_avg!2  = "lat"
Qh_avg!3  = "lon"
Qh_avg&lat= CABLE_lat
Qh_avg&lon= CABLE_lon

deltaT_avg!0  = "case"
deltaT_avg!1  = "time_period"
deltaT_avg!2  = "lat"
deltaT_avg!3  = "lon"
deltaT_avg&lat= CABLE_lat
deltaT_avg&lon= CABLE_lon

Fwsoil_avg!0  = "case"
Fwsoil_avg!1  = "time_period"
Fwsoil_avg!2  = "lat"
Fwsoil_avg!3  = "lon"
Fwsoil_avg&lat= CABLE_lat
Fwsoil_avg&lon= CABLE_lon

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
end if

do case_num = 0, case_sum -1
   day_s = 0
   do year = year_s, year_e
      print("year = " +year)
      if ( mod(year,4) .eq. 0) then
        doy = 366
      else
        doy = 365
      end if

      filename = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
      f        = addfile (filename,"r")

      Qle(day_s:day_s+doy-1,:,:)    = f->Qle(:,{A:B},{C:D})
      Qh(day_s:day_s+doy-1,:,:)     = f->Qh(:,{A:B},{C:D})
      Fwsoil(day_s:day_s+doy-1,:,:) = f->Fwsoil(:,{A:B},{C:D})
      deltaT(day_s:day_s+doy-1,:,:) = f->VegT(:,{A:B},{C:D}) - f->Tair(:,{A:B},{C:D})
      day_s = day_s + doy
      delete(f)
   end do ; year

   if heatwave .eq. "hw-90th" then
      Qle    = where(HW .eq. 1, Qle, default_val)
      Qh     = where(HW .eq. 1, Qh, default_val)
      Fwsoil = where(HW .eq. 1, Fwsoil, default_val)
      deltaT = where(HW .eq. 1, deltaT, default_val)
   end if

   do count = 0, 3
     Qle_avg(case_num, count,:,:)    = dim_avg_n(Qle(ts_s(count):ts_e(count),:,:), 0)
     Qh_avg(case_num, count,:,:)     = dim_avg_n(Qh(ts_s(count):ts_e(count),:,:), 0)
     Fwsoil_avg(case_num, count,:,:) = dim_avg_n(Fwsoil(ts_s(count):ts_e(count),:,:), 0)
     deltaT_avg(case_num, count,:,:) = dim_avg_n(deltaT(ts_s(count):ts_e(count),:,:), 0)
   end do
end do ; case_num

delete(Qle)
delete(Qh)
delete(Fwsoil)
delete(deltaT)
delete(HW)

; ; ================================ READ LST ==================================
; ; use 2000-2019
; LST_ET = path_LST+"/E_1980_2018_LST_v3.3a_YR.nc"
; G_ET     = addfile (LST_ET,"r")
;
; E_tmp    = G_ET->E(:,{C:D},{A:B})
; E_tmp    = doubletofloat(E_tmp)
;
; Evap_LST = E_tmp(time|20:38,lat|:, lon|:)
;
; delete(G_ET)

; ================================= Plotting ===================================
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")

res                         = True
res@gsnDraw                 = False  ; Don't draw plots
res@gsnFrame                = False  ; ����ҳ
res@gsnMaximize             = False

res@vpHeightF 	            = 0.4               ; change aspect ratio of plot
res@vpWidthF 	              = 0.8
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

res@cnFillMode             = "RasterFill"            ; Raster Mode
res@cnFillOn               = True                            ;��ɫ
res@cnLinesOn              = False
res@cnLineLabelsOn         = False
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode   = "ExplicitLevels"

res@lbLabelBarOn           = True
res@lbBoxLinesOn           = True                       ;�ر�lbar box ����
; res@lbOrientation          = "vertical"
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

pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False
pres@cnLevelSelectionMode = "ExplicitLevels"
pres@lbBoxLinesOn       = True
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01
pres@gsnPanelYWhiteSpacePercent = 0.5
pres@gsnPanelXWhiteSpacePercent = 0.5

txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.03

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

Qle_CTL_lvl    = (/10., 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140/)
Fw_CTL_lvl     = fspan(0.1,0.9,9)
Qh_CTL_lvl     = (/20., 40, 60, 80, 100, 120, 140, 160, 180, 200, 220/)
deltaT_CTL_lvl = fspan(0.5,5,10)

; diff_color = "BlAqGrWh2YeOrReVi22"; "ViBlGrWhYeOrRe";
Qle_diff_lvl    = (/-40.,-35.,-30.,-25.,-20., -15., -10., -5., -2.5, 0, 2.5, 5., 10., 15., 20.,25.,30.,35.,40./)
Fw_diff_lvl     = fspan(-0.4,0.4,17)
Qh_diff_lvl     = Qle_diff_lvl
deltaT_diff_lvl = (/-1., -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1./)

plots = new(case_sum*4,graphic)

do cnt = 1,2
  do case_num = 0,case_sum-1
    if case_num .eq. 0 then
      res@lbLabelBarOn  = True

      ; deltaT
      res@cnLevels      = deltaT_CTL_lvl
      res@cnFillPalette = "precip3_16lev";"WhiteBlueGreenYellowRed";"vegetation_ClarkU";"NCV_blue_red";"NOC_ndvi";"MPL_RdYlGn";"WhiteBlueGreenYellowRed"
      plots(case_num*4) = gsn_csm_contour_map(wks,deltaT_avg(0,cnt,:,:),res)
      delete(res@cnLevels)

      ; Qh
      res@cnLevels      = Qh_CTL_lvl
      plots(case_num*4+1) = gsn_csm_contour_map(wks,Qh_avg(0,cnt,:,:),res)
      delete(res@cnLevels)
      delete(res@cnFillPalette)

      ; Qle
      res@cnLevels      = Qle_CTL_lvl
      ;cmap = read_colormap_file("WhiteBlueGreenYellowRed");("NEO_div_vegetation_a")("ViBlGrWhYeOrRe");
      res@cnFillPalette = "spread_15lev";cmap(::-1,:)
      plots(case_num*4+2) = gsn_csm_contour_map(wks,Qle_avg(0,cnt,:,:),res)
      delete(res@cnLevels)

      ; Fwsoil
      res@cnLevels      = Fw_CTL_lvl
      plots(case_num*4+3) = gsn_csm_contour_map(wks,Fwsoil_avg(0,cnt,:,:),res)
      delete(res@cnLevels)

    else if case_num .lt. 3 then
      res@lbLabelBarOn  = False

      delete(res@cnFillPalette)
      res@cnLevels      = deltaT_diff_lvl
      res@cnFillPalette = "BlAqGrWh2YeOrReVi22";"BlAqGrWh2YeOrReVi22"
      tmp  = deltaT_avg(case_num,cnt,:,:)
      tmp  = tmp - deltaT_avg(0,cnt,:,:)
      plots(case_num*4) = gsn_csm_contour_map(wks,tmp,res)
      delete(tmp)
      delete(res@cnLevels)

      res@cnLevels      = Qh_diff_lvl
      tmp  = Qh_avg(case_num,cnt,:,:)
      tmp  = tmp - Qh_avg(0,cnt,:,:)
      plots(case_num*4+1) = gsn_csm_contour_map(wks,tmp,res)
      delete(tmp)
      delete(res@cnLevels)

      ; delete(cmap)
      delete(res@cnFillPalette)
      cmap = read_colormap_file("BlAqGrWh2YeOrReVi22");("NEO_div_vegetation_a")("ViBlGrWhYeOrRe");
      res@cnFillPalette = cmap(::-1,:)

      res@cnLevels      = Qle_diff_lvl
      tmp  = Qle_avg(case_num,cnt,:,:)
      tmp  = tmp - Qle_avg(0,cnt,:,:)
      plots(case_num*4+2) = gsn_csm_contour_map(wks,tmp,res)
      delete(tmp)
      delete(res@cnLevels)

      res@cnLevels      = Fw_diff_lvl
      tmp  = Fwsoil_avg(case_num,cnt,:,:)
      tmp  = tmp - Fwsoil_avg(0,cnt,:,:)
      plots(case_num*4+3) = gsn_csm_contour_map(wks,tmp,res)
      delete(tmp)
      delete(res@cnLevels)

    else if case_num .eq. 3 then
      res@lbLabelBarOn  = True

      delete(res@cnFillPalette)
      res@cnLevels      = deltaT_diff_lvl
      res@cnFillPalette = "BlAqGrWh2YeOrReVi22";"BlAqGrWh2YeOrReVi22"
      tmp  = deltaT_avg(case_num,cnt,:,:)
      tmp  = tmp - deltaT_avg(1,cnt,:,:)
      plots(case_num*4) = gsn_csm_contour_map(wks,tmp,res)
      delete(tmp)
      delete(res@cnLevels)

      res@cnLevels      = Qh_diff_lvl
      tmp  = Qh_avg(case_num,cnt,:,:)
      tmp  = tmp - Qh_avg(1,cnt,:,:)
      plots(case_num*4+1) = gsn_csm_contour_map(wks,tmp,res)
      delete(tmp)
      delete(res@cnLevels)

      ; delete(cmap)
      delete(res@cnFillPalette)
      cmap = read_colormap_file("BlAqGrWh2YeOrReVi22");("NEO_div_vegetation_a")("ViBlGrWhYeOrRe");
      res@cnFillPalette = cmap(::-1,:)

      res@cnLevels      = Qle_diff_lvl
      tmp  = Qle_avg(case_num,cnt,:,:)
      tmp  = tmp - Qle_avg(1,cnt,:,:)
      plots(case_num*4+2) = gsn_csm_contour_map(wks,tmp,res)
      delete(tmp)
      delete(res@cnLevels)

      res@cnLevels      = Fw_diff_lvl
      tmp  = Fwsoil_avg(case_num,cnt,:,:)
      tmp  = tmp - Fwsoil_avg(1,cnt,:,:)
      plots(case_num*4+3) = gsn_csm_contour_map(wks,tmp,res)
      delete(tmp)
      delete(res@cnLevels)
    end if
    end if
    end if
  end do

  do count = 0, 15
    txid1   = gsn_create_text(wks, order_name(count), txres)
    annoid1 = gsn_add_annotation(plots(count), txid1, amres)
  end do

  gsn_panel(wks,(/plots/),(/4,4/),pres)
  delete(res@cnFillPalette)
end do ; cnt
exit
end
