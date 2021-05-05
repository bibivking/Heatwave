
; MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
; Author  : Mengyuan Mu
; Version : 1.0 (22.03.2021)"
; Email   : mu.mengyuan815@gmail.com
; WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

;*******************************************
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

; =================================== Path =====================================
case_name   = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr"/)

path_mum    = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
case_sum    = dimsizes(case_name)
path        = new(case_sum,"string")
pic_message = "CTL_vs_FD"
order_name  = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"/) ;
categorize  = True
colorbar    = "MPL_BrBG";"CBR_drywet";"BlueDarkOrange18";"MPL_BrBG"; "NEO_div_vegetation_b";"NEO_div_vegetation_c";"MPL_BrBG" "CBR_wet"
cmap_tmp    = read_colormap_file(colorbar)
cmap        = new((/19, 4/),float)
cmap(0,:)   = cmap_tmp(0,:)
cmap(1,:)   = cmap_tmp(7,:)
cmap(2,:)   = cmap_tmp(15,:)
cmap(3,:)   = cmap_tmp(21,:)
cmap(4,:)   = cmap_tmp(27,:)
cmap(5,:)   = cmap_tmp(35,:)
cmap(6,:)   = cmap_tmp(37,:)
cmap(7,:)   = cmap_tmp(44,:)
cmap(8,:)   = cmap_tmp(55,:)
cmap(9,:)   = cmap_tmp(63,:)
cmap(10,:)   = cmap_tmp(72,:)
cmap(11,:)   = cmap_tmp(79,:)
cmap(12,:)   = cmap_tmp(86,:)
cmap(13,:)   = cmap_tmp(95,:)
cmap(14,:)   = cmap_tmp(101,:)
cmap(15,:)   = cmap_tmp(111,:)
cmap(16,:)   = cmap_tmp(116,:)
cmap(17,:)   = cmap_tmp(122,:)
cmap(18,:)   = cmap_tmp(127,:)

do case_num = 0, case_sum -1
   path(case_num) = path_mum+case_name(case_num)
end do

; ================================= Predefine ==================================
year_s   = 2017
year_e   = 2019
year_sum = year_e-year_s+1
row_sum  = 2
col_sum  = 3

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

; ================================== OPTIONS ===================================
plot_type = "all"
scale     = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"

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

filename  = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f         = addfile(filename,"r")

lat_sum   = dimsizes(f->Rainf(0,{A:B},0))
lon_sum   = dimsizes(f->Rainf(0,0,{C:D}))
lat_size  = dimsizes(f->y)
lon_size  = dimsizes(f->x)
lat       = f->y({A:B})
lon       = f->x({C:D})
default_val = -1.e+33
delete(f)

Evap      = new((/case_sum, year_sum, lat_sum, lon_sum/),float, default_val)
Trans     = new((/case_sum, year_sum, lat_sum, lon_sum/),float, default_val)
Esoil     = new((/case_sum, year_sum, lat_sum, lon_sum/),float, default_val)
Fwsoil    = new((/case_sum, year_sum, lat_sum, lon_sum/),float, default_val)
Q2botsoil = new((/case_sum, year_sum, lat_sum, lon_sum/),float, default_val) ; Qrec from the aqufier to bottom soil layer
Q2aquifer = new((/case_sum, year_sum, lat_sum, lon_sum/),float, default_val) ; Qrec from bottom soil layer to aquifer

Evap      = 0.
Trans     = 0.
Esoil     = 0.
Fwsoil    = 0.
Q2botsoil = 0.
Q2aquifer = 0.

do case_num = 0, case_sum -1
   do year = year_s, year_e

      print("year = " +year)

      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if

      filename = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
      f        = addfile (filename,"r")

      Evap(case_num,year-year_s,:,:)     = dim_sum_n_Wrap(f->Evap(:,{A:B},{C:D}),0)
      Trans(case_num,year-year_s,:,:)    = dim_sum_n_Wrap(f->TVeg(:,{A:B},{C:D}),0)
      Esoil(case_num,year-year_s,:,:)    = dim_sum_n_Wrap(f->ESoil(:,{A:B},{C:D}),0)
      Fwsoil(case_num,year-year_s,:,:)   = dim_avg_n_Wrap(f->Fwsoil(:,{A:B},{C:D}),0)

      Qrec_tmp = f->Qrecharge(:,{A:B},{C:D})
      Q2botsoil(case_num,year-year_s,:,:) = dim_sum_n_Wrap(where(Qrec_tmp .le. 0, Qrec_tmp, 0),0)
      Q2aquifer(case_num,year-year_s,:,:) = dim_sum_n_Wrap(where(Qrec_tmp .ge. 0, Qrec_tmp, 0),0)

      delete(f)
      delete(Qrec_tmp)
   end do ; year
end do ; case_num

Evap     = Evap*3600.*24.
Trans    = Trans*3600.*24.
Esoil    = Esoil*3600.*24.
Q2botsoil= Q2botsoil*3600.*24.*(-1.)
Q2aquifer= Q2aquifer*3600.*24.

Evap!0           = "case"
Evap!1           = "year"
Evap!2           = "lat"
Evap!3           = "lon"
Evap&lat         = lat
Evap&lon         = lon

Trans!0           = "case"
Trans!1           = "year"
Trans!2           = "lat"
Trans!3           = "lon"
Trans&lat         = lat
Trans&lon         = lon

Esoil!0           = "case"
Esoil!1           = "year"
Esoil!2           = "lat"
Esoil!3           = "lon"
Esoil&lat         = lat
Esoil&lon         = lon

Q2botsoil!0      = "case"
Q2botsoil!1      = "year"
Q2botsoil!2      = "lat"
Q2botsoil!3      = "lon"
Q2botsoil&lat    = lat
Q2botsoil&lon    = lon

Q2aquifer!0      = "case"
Q2aquifer!1      = "year"
Q2aquifer!2      = "lat"
Q2aquifer!3      = "lon"
Q2aquifer&lat    = lat
Q2aquifer&lon    = lon

; ============================== Max Fwsoil_diff ===============================
delete(case_name)
case_name  = (/"CTL-30x3+30yr","FREE_DRAIN"/)
HW_day_s   = (/ 15, 22/)
HW_day_e   = (/ 19, 27/)
event_sum  = dimsizes(HW_day_s)
case_sum   = dimsizes(case_name)

day_tot    = 0
do event_num = 0, event_sum-1
  day_tot = day_tot + (HW_day_e(event_num)-HW_day_s(event_num)+1)
end do

path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = path_mum+case_name(case_num)+"/outputs-subday"
end do

Fwsoil_subday   = new((/case_sum,lat_sum,lon_sum/),float,default_val)
Fwsoil_subday   = 0.

Fwsoil_diff_max = new((/case_sum,lat_sum,lon_sum/),float,default_val)
Fwsoil_diff_max = 0

do case_num = 0, case_sum -1
  print("case_num = " +case_num)
  filename1 = path(case_num)+"/cable_out_2019_SE_Aus.nc"
  f1 = addfile (filename1,"r")

  Rain_tmp     = f1->Rainf(0:60*8,{A:B},{C:D})
  Rain_tmp     = Rain_tmp*3*3600
  Fwsoil_tmp   = f1->Fwsoil(0:60*8,{A:B},{C:D})
  delete(f1)
  delete(filename1)

  do event_num = 0,event_sum-1
    do day = HW_day_s(event_num),HW_day_e(event_num)
      ts = (day-1)*8+7 ;6-6pm ;7 - 9pm
      Rain = dim_sum_n_Wrap(Rain_tmp(ts:ts+3,:,:),0)
      ; 3am - 9pm
      Fwsoil_subday(case_num,:,:)   = where(Rain .lt. 0.01, \
                                            Fwsoil_tmp(ts+2,:,:)-Fwsoil_tmp(ts,:,:), default_val)
      Fwsoil_diff_max(case_num,:,:) = where(Fwsoil_subday(case_num,:,:) .gt. Fwsoil_diff_max(case_num,:,:),\
                                            Fwsoil_subday(case_num,:,:), Fwsoil_diff_max(case_num,:,:))
      delete(Rain)
    end do
  end do
end do ; case_num

Fwsoil_diff_max!0      = "time"
Fwsoil_diff_max!1      = "lat"
Fwsoil_diff_max!2      = "lon"
Fwsoil_diff_max&lat    = lat
Fwsoil_diff_max&lon    = lon

; ================================ Read WTD ====================================
WTD         = new((/241,281/),float)
WTD         = 0.

do year  = year_s,year_e
  filename_hw = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/CTL-30x3+30yr/outputs-daily/cable_out_"+year+"_SE_Aus.nc"
  f           = addfile (filename_hw,"r")
  WTD         = WTD + dim_avg_n_Wrap(f->WatTable(:,{A:B},{C:D}),0)
end do

WTD = WTD/year_sum
WTD@_FillValue = f->WatTable@_FillValue
WTD!0 = "lat"
WTD!1 = "lon"
WTD&lat = f->y ;fspan(-40,-28,241)
WTD&lon = f->x ;fspan(140,154,281)
printVarSummary(WTD)

; ================================= Plotting ===================================
res                         = True
res@gsnDraw                 = True  ; Don't draw plots
res@gsnFrame                = True  ; ����ҳ
res@gsnMaximize             = True

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

res@lbLabelBarOn            = True
res@lbBoxLinesOn            = True                       ;�ر�lbar box ����
res@lbTitleFont             = 25
res@lbLabelFont             = 25
res@lbTitleFontHeightF      = 0.015
res@lbLabelFontHeightF      = 0.015
res@lbLabelAngleF           = 90
; res@lbOrientation           = "Vertical"
res@lbTitlePosition         = "Bottom"
res@pmLabelBarHeightF       = 0.1
res@pmLabelBarParallelPosF  = 0.5
res@pmLabelBarOrthogonalPosF= 0.12

res@tmXBLabelFont           = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont           = 25
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF    = 0.015
res@tmYLLabelFontHeightF    = 0.015
; res@mpDataBaseVersion       = "HighRes"
res@tiMainFontHeightF       = 0.02

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
  res@tmYLLabels      = (/"40~S~o~N~N","35~S~o~N~N","30~S~o~N~","25~S~o~N~S","20~S~o~N~N","15~S~o~N~N","10~S~o~N~"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/115,120,125,130,135,140,145,150/)
  res@tmXBLabels      = (/"115~S~o~N~W","120~S~o~N~W","125~S~o~N~","130~S~o~N~W","135~S~o~N~W","140~S~o~N~","145~S~o~N~W","150~S~o~N~"/)
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

txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.05

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

; ________________ add region average ________________
txres1                  = True
txres1@txPerimOn        = False
txres1@txFont           = 25
txres1@txFontHeightF    = 0.023
amres1                  = True
amres1@amParallelPosF   = 0.45    ; This is the right edge of the plot.
amres1@amOrthogonalPosF = 0.45    ; This is the bottom edge of the plot.
amres1@amJust           = "BottomRight"

; =========================== plot Evaptranspiration ===========================
pic = "./plots/Fig2_spitial_map_E_Qrec_Fwmax_"+pic_message

; higher resolution
wks_type          = "pdf"
wks_type@wkWidth  = 2500
wks_type@wkHeight = 2500
wks               = gsn_open_wks(wks_type,pic)
gsn_define_colormap(wks,"rainbow") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

plots             = new(6,graphic)

res@cnFillOn      = True
res@cnMissingValFillColor = "gray80"
res@cnFillBackgroundColor = "gray80"

res@cnLevels      = (/\
                      -400,-300,-200,-150,-100,-50,-25,-5, \
                      5,25,50,100,150,200,300,400         \
                      /)

res@lbLabelStrings =  (/ "-400","-300","-200","-150","-100","-50","-25","-5", \
                         "5","25","50","100","150","200","300","400"        \
                       /)


; -------------------------------------------
if categorize then
  tmp_WTD_class = new((/6, lat_sum, lon_sum/),float, default_val)
  tmp_WTD_class@_FillValue = -9999.
  tmp_WTD_class!0 = "WTD_lvl"
  tmp_WTD_class!1 = "lat"
  tmp_WTD_class!2 = "lon"
  tmp_WTD_class&lat = f->y ;fspan(-40,-28,241)
  tmp_WTD_class&lon = f->x ;fspan(140,154,281)
end if
; -------------------------------------------

res@cnFillPalette = cmap
; res@lbTitleString = "~F5~D"+"~F25~E (mm yr~S~-1~N~)"
res@tiMainString  = "~F5~D"+"~F25~E (mm yr~S~-1~N~)"
tmp               = dim_avg_n_Wrap(Evap(0,:,:,:),0)
tmp               = tmp - dim_avg_n_Wrap(Evap(1,:,:,:),0)
tmp@_FillValue    = -9999.

tmp_WTD_0         = avg(where(WTD .ge. 0.0 .and. WTD .lt. 5.0, tmp, tmp@_FillValue))
tmp_WTD_5         = avg(where(WTD .ge. 5.0 .and. WTD .lt. 10.0, tmp, tmp@_FillValue))
tmp_WTD_10        = avg(where(WTD .ge. 10.0, tmp, tmp@_FillValue))
context           = "~C~ WTD  "+"~F5~D"+"~F25~E"               \
                  + "~C~  0-5m "+sprintf("%4.0f",tmp_WTD_0)  \
                  + "~C~5-10m "+sprintf("%4.0f",tmp_WTD_5)  \
                  + "~C~ >10m "+sprintf("%4.0f",tmp_WTD_10)

plots(0)          = gsn_csm_contour_map(wks,tmp,res)
txid2             = gsn_create_text(wks, context, txres1)
annoid2           = gsn_add_annotation(plots(0), txid2, amres1)

; ------------------------------------------
if categorize then
  tmp_WTD_class(0,:,:) = where(tmp .lt. 25,tmp,tmp@_FillValue)
  tmp_WTD_class(1,:,:) = where(tmp .ge. 25 .and. tmp .lt. 50,tmp,tmp@_FillValue)
  tmp_WTD_class(2,:,:) = where(tmp .ge. 50 .and. tmp .lt. 100,tmp,tmp@_FillValue)
  tmp_WTD_class(3,:,:) = where(tmp .ge. 100 .and. tmp .lt. 200,tmp,tmp@_FillValue)
  tmp_WTD_class(4,:,:) = where(tmp .ge. 200 .and. tmp .lt. 500,tmp,tmp@_FillValue)
  tmp_WTD_class(5,:,:) = where(tmp .ge. 500,tmp,tmp@_FillValue)

  res@gsnLeftString = "E < 25"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(0,:,:), res)
  res@gsnLeftString = "25 < E < 50"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(1,:,:), res)
  res@gsnLeftString = "50 < E < 100"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(2,:,:), res)
  res@gsnLeftString = "100 < E < 200"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(3,:,:), res)
  res@gsnLeftString = "200 < E < 500"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(4,:,:), res)
  res@gsnLeftString = "500 < E"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(5,:,:), res)
  res@gsnLeftString =""
end if
; ------------------------------------------

delete(tmp)
delete(tmp_WTD_0)
delete(tmp_WTD_5)
delete(tmp_WTD_10)
delete(context)

; res@lbTitleString = "~F5~D"+"~F25~Et (mm yr~S~-1~N~)"
res@tiMainString = "~F5~D"+"~F25~Et (mm yr~S~-1~N~)"
tmp               = dim_avg_n_Wrap(Trans(0,:,:,:),0)
tmp               = tmp - dim_avg_n_Wrap(Trans(1,:,:,:),0)
tmp@_FillValue    = -9999.

tmp_WTD_0         = avg(where(WTD .ge. 0.0 .and. WTD .lt. 5.0, tmp, tmp@_FillValue))
tmp_WTD_5         = avg(where(WTD .ge. 5.0 .and. WTD .lt. 10.0, tmp, tmp@_FillValue))
tmp_WTD_10        = avg(where(WTD .ge. 10.0, tmp, tmp@_FillValue))
context           = "~C~ WTD  "+"~F5~D"+"~F25~Et"               \
                  + "~C~  0-5m "+sprintf("%4.0f",tmp_WTD_0)  \
                  + "~C~5-10m "+sprintf("%4.0f",tmp_WTD_5)  \
                  + "~C~ >10m "+sprintf("%4.0f",tmp_WTD_10)

plots(1)          = gsn_csm_contour_map(wks,tmp,res)
txid2             = gsn_create_text(wks, context, txres1)
annoid2           = gsn_add_annotation(plots(1), txid2, amres1)

; ------------------------------------------
  if categorize then
  tmp_WTD_class(0,:,:) = where(tmp .lt. 25,tmp,tmp@_FillValue)
  tmp_WTD_class(1,:,:) = where(tmp .ge. 25 .and. tmp .lt. 50,tmp,tmp@_FillValue)
  tmp_WTD_class(2,:,:) = where(tmp .ge. 50 .and. tmp .lt. 100,tmp,tmp@_FillValue)
  tmp_WTD_class(3,:,:) = where(tmp .ge. 100 .and. tmp .lt. 200,tmp,tmp@_FillValue)
  tmp_WTD_class(4,:,:) = where(tmp .ge. 200 .and. tmp .lt. 500,tmp,tmp@_FillValue)
  tmp_WTD_class(5,:,:) = where(tmp .ge. 500,tmp,tmp@_FillValue)

  res@gsnLeftString = "Et < 25"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(0,:,:), res)
  res@gsnLeftString = "25 < Et < 50"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(1,:,:), res)
  res@gsnLeftString = "50 < Et < 100"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(2,:,:), res)
  res@gsnLeftString = "100 < Et < 200"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(3,:,:), res)
  res@gsnLeftString = "200 < Et < 500"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(4,:,:), res)
  res@gsnLeftString = "500 < Et"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(5,:,:), res)
  res@gsnLeftString =""
end if
; ------------------------------------------

delete(tmp)
delete(tmp_WTD_0)
delete(tmp_WTD_5)
delete(tmp_WTD_10)
delete(context)

; res@lbTitleString = "~F5~D"+"~F25~Es (mm yr~S~-1~N~)"
res@tiMainString = "~F5~D"+"~F25~Es (mm yr~S~-1~N~)"
tmp               = dim_avg_n_Wrap(Esoil(0,:,:,:),0)
tmp               = tmp - dim_avg_n_Wrap(Esoil(1,:,:,:),0)
tmp@_FillValue    = -9999.

tmp_WTD_0         = avg(where(WTD .ge. 0.0 .and. WTD .lt. 5.0, tmp, tmp@_FillValue))
tmp_WTD_5         = avg(where(WTD .ge. 5.0 .and. WTD .lt. 10.0, tmp, tmp@_FillValue))
tmp_WTD_10        = avg(where(WTD .ge. 10.0, tmp, tmp@_FillValue))
context           = "~C~ WTD  "+"~F5~D"+"~F25~Es"               \
                  + "~C~  0-5m "+sprintf("%4.0f",tmp_WTD_0)  \
                  + "~C~5-10m "+sprintf("%4.0f",tmp_WTD_5)  \
                  + "~C~ >10m "+sprintf("%4.0f",tmp_WTD_10)

amres1@amParallelPosF   = 0.485    ; This is the right edge of the plot.
amres1@amOrthogonalPosF = 0.45    ; This is the bottom edge of the plot.
plots(2)          = gsn_csm_contour_map(wks,tmp,res)
txid2             = gsn_create_text(wks, context, txres1)
annoid2           = gsn_add_annotation(plots(2), txid2, amres1)

; ------------------------------------------
if categorize then
  tmp_WTD_class(0,:,:) = where(tmp .lt. 5 .and. tmp .lt. 25, tmp,tmp@_FillValue)
  tmp_WTD_class(1,:,:) = where(tmp .ge. 25 .and. tmp .lt. 50,tmp,tmp@_FillValue)
  tmp_WTD_class(2,:,:) = where(tmp .ge. 50 .and. tmp .lt. 100,tmp,tmp@_FillValue)
  tmp_WTD_class(3,:,:) = where(tmp .ge. 100 .and. tmp .lt. 200,tmp,tmp@_FillValue)
  tmp_WTD_class(4,:,:) = where(tmp .ge. 200 .and. tmp .lt. 500,tmp,tmp@_FillValue)
  tmp_WTD_class(5,:,:) = where(tmp .ge. 500,tmp,tmp@_FillValue)

  res@gsnLeftString = "5 < Es < 25"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(0,:,:), res)
  res@gsnLeftString = "25 < Es < 50"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(1,:,:), res)
  res@gsnLeftString = "50 < Es < 100"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(2,:,:), res)
  res@gsnLeftString = "100 < Es < 200"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(3,:,:), res)
  res@gsnLeftString = "200 < Es < 500"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(4,:,:), res)
  res@gsnLeftString = "500 < Es"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(5,:,:), res)
  res@gsnLeftString =""
end if
; ------------------------------------------

delete(res@cnFillPalette)
delete(tmp)
delete(tmp_WTD_0)
delete(tmp_WTD_5)
delete(tmp_WTD_10)
delete(context)


; =============================== Plot Dr Qrec ================================
; Q2aquifer
print("plot Q2aquifer")

res@cnFillPalette = cmap(::-1,:)
; res@lbTitleString = "~F5~D"+"~F25~Dr (mm yr~S~-1~N~)"
res@tiMainString = "~F5~D"+"~F25~Dr (mm yr~S~-1~N~)"
tmp               = dim_avg_n_Wrap(Q2aquifer(0,:,:,:),0)
tmp               = tmp - dim_avg_n_Wrap(Q2aquifer(1,:,:,:),0)
tmp@_FillValue    = -9999.

tmp_WTD_0         = avg(where(WTD .ge. 0.0 .and. WTD .lt. 5.0, tmp, tmp@_FillValue))
tmp_WTD_5         = avg(where(WTD .ge. 5.0 .and. WTD .lt. 10.0, tmp, tmp@_FillValue))
tmp_WTD_10        = avg(where(WTD .ge. 10.0, tmp, tmp@_FillValue))
context           = "~C~ WTD  "+"~F5~D"+"~F25~Dr"           \
                  + "~C~  0-5m "+sprintf("%4.0f",tmp_WTD_0)  \
                  + "~C~5-10m "+sprintf("%4.0f",tmp_WTD_5)  \
                  + "~C~ >10m "+sprintf("%4.0f",tmp_WTD_10)

plots(3)          = gsn_csm_contour_map(wks,tmp,res)
txid2             = gsn_create_text(wks, context, txres1)
annoid2           = gsn_add_annotation(plots(3), txid2, amres1)

; ------------------------------------------
if categorize then
  tmp_WTD_class(0,:,:) = where(tmp .ge. 50,tmp,tmp@_FillValue)
  tmp_WTD_class(1,:,:) = where(tmp .ge. 5 .and. tmp .lt. 50, tmp,tmp@_FillValue)
  tmp_WTD_class(2,:,:) = where(tmp .ge. -50 .and. tmp .lt. -5,   tmp,tmp@_FillValue)
  tmp_WTD_class(3,:,:) = where(tmp .ge. -200 .and. tmp .lt. -50,  tmp,tmp@_FillValue)
  tmp_WTD_class(4,:,:) = where(tmp .ge. -500 .and. tmp .lt. -200,tmp,tmp@_FillValue)
  tmp_WTD_class(5,:,:) = where(tmp .lt. -500,tmp,tmp@_FillValue)

  res@gsnLeftString = "50 < Dr"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(0,:,:), res)
  res@gsnLeftString = "5 < Dr < 50"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(1,:,:), res)
  res@gsnLeftString = "-50 < Dr < -5"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(2,:,:), res)
  res@gsnLeftString = "-200 < Dr < -50"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(3,:,:), res)
  res@gsnLeftString = "-500 < Dr < -200"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(4,:,:), res)
  res@gsnLeftString = "Dr < -500"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(5,:,:), res)

  ; res@gsnLeftString = "Dr < -200"
  ; tmp_WTD_class(5,:,:) = where(tmp .lt. -200, tmp,tmp@_FillValue)
  ; plot              = gsn_csm_contour_map(wks,tmp_WTD_class(5,:,:), res)

  res@gsnLeftString =""
end if
; ------------------------------------------

delete(res@cnFillPalette)
delete(tmp)
delete(tmp_WTD_0)
delete(tmp_WTD_5)
delete(tmp_WTD_10)
delete(context)

; Q2botsoil
print("plot Q2botsoil")

; res@cnLevels      = (/\
;                       -600,-500,-400,-300,-200,-100,-50,-40,-30,-20,-10, \
;                       0,10,20,30,40,50,100,200,300,400,500,600          \
;                       /)

; res@lbTitleString = "~F5~D"+"~F25~Qrec (mm yr~S~-1~N~)"
res@tiMainString = "~F5~D"+"~F25~Qrec (mm yr~S~-1~N~)"
res@cnFillPalette = cmap
tmp               = dim_avg_n_Wrap(Q2botsoil(0,:,:,:),0)
tmp               = tmp - dim_avg_n_Wrap(Q2botsoil(1,:,:,:),0)
tmp@_FillValue    = -9999.

tmp_WTD_0         = avg(where(WTD .ge. 0.0 .and. WTD .lt. 5.0, tmp, tmp@_FillValue))
tmp_WTD_5         = avg(where(WTD .ge. 5.0 .and. WTD .lt. 10.0, tmp, tmp@_FillValue))
tmp_WTD_10        = avg(where(WTD .ge. 10.0, tmp, tmp@_FillValue))
context           = "~C~ WTD  "+"~F5~D"+"~F25~Qrec"               \
                  + "~C~  0-5m "+sprintf("%4.0f",tmp_WTD_0)  \
                  + "~C~5-10m "+sprintf("%4.0f",tmp_WTD_5)  \
                  + "~C~ >10m "+sprintf("%4.0f",tmp_WTD_10)

plots(4)          = gsn_csm_contour_map(wks,tmp,res)
txid2             = gsn_create_text(wks, context, txres1)
annoid2           = gsn_add_annotation(plots(4), txid2, amres1)

; ------------------------------------------
if categorize then
  tmp_WTD_class(0,:,:) = where(tmp .lt. 0,tmp,tmp@_FillValue)
  tmp_WTD_class(1,:,:) = where(tmp .ge. 0 .and. tmp .lt. 25,tmp,tmp@_FillValue)
  tmp_WTD_class(2,:,:) = where(tmp .ge. 25 .and. tmp .lt. 50,tmp,tmp@_FillValue)
  tmp_WTD_class(3,:,:) = where(tmp .ge. 50 .and. tmp .lt. 100,tmp,tmp@_FillValue)
  tmp_WTD_class(4,:,:) = where(tmp .ge. 100 .and. tmp .lt. 200,tmp,tmp@_FillValue)
  tmp_WTD_class(5,:,:) = where(tmp .ge. 200,tmp,tmp@_FillValue)

  res@gsnLeftString = "Qrec < 0"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(0,:,:), res)
  res@gsnLeftString = "0 < Qrec < 25"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(1,:,:), res)
  res@gsnLeftString = "25 < Qrec < 50"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(2,:,:), res)
  res@gsnLeftString = "50 < Qrec < 100"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(3,:,:), res)
  res@gsnLeftString = "100 < Qrec < 200"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(4,:,:), res)
  res@gsnLeftString = "200 < Qrec"
  plot              = gsn_csm_contour_map(wks,tmp_WTD_class(5,:,:), res)
  res@gsnLeftString =""
end if
; ------------------------------------------

delete(res@cnFillPalette)
delete(res@cnLevels)
delete(tmp)
delete(tmp_WTD_0)
delete(tmp_WTD_5)
delete(tmp_WTD_10)
delete(context)


; ============================ Plot Fwsoil_diff_max ============================
delete(res@lbLabelStrings)
; Fwsoil
print("plot Fwsoil_diff_max")
; res@lbTitleString = "~F5~Db"
res@tiMainString = "~F5~Db ~F25~(-)"
res@cnFillPalette = cmap(9:,:)
res@cnLevels      = fspan(0.005,0.045,9)
plots(5) = gsn_csm_contour_map(wks,Fwsoil_diff_max(0,:,:),res)

do counter = 0, 5
  txid1   = gsn_create_text(wks, order_name(counter), txres)
  annoid1 = gsn_add_annotation(plots(counter), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/2,3/),pres)
end
