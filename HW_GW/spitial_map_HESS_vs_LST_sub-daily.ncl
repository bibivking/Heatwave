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
PFT        = "all" ; "tree"
scale      = "SE-AU"
year_s     = 2019
year_e     = 2019
diff       = False

HW_day_s   = (/ 1 /) ;15, 22/)
HW_day_e   = (/ 31/) ;19, 27/)
event_sum  = dimsizes(HW_day_s)

day_tot    = 0
do event_num = 0, event_sum-1
  day_tot = day_tot + (HW_day_e(event_num)-HW_day_s(event_num)+1)
end do

if diff .eq. True then
  message    = "LST_subday_diff"
else
  message    = "LST_subday"
end if
pic_name   = "HESS_"+message+"_"+scale+"_"+year_s+"-"+year_e

path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
case_name = (/"CTL-30x3+30yr","FREE_DRAIN","CTL-30x3+30yr_2019-DeepRoot"/)
order_name  = (/"(a)","(b)","(c)","(d)",\
                "(e)","(f)","(g)","(h)",\
                "(i)","(j)","(k)","(l)",\
                "(m)","(n)","(o)","(p)",\
                "(q)","(r)","(s)","(t)"/)

; 0, 1, 2, 3 , 4 , 5,  6,  7
; 0, 3, 6, 9, 12, 15, 18, 21, 24

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
f1  = addfile (filename1,"r")
lat = f1->y({A:B})
lon = f1->x({C:D})

PFT_mask = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
default_val = f1->Qle@_FillValue

delete(filename1)
delete(f1)

VegT     = new((/case_sum,2,31,dimsizes(lat),dimsizes(lon)/),float,default_val)
VegT     = 0.

; ================================== VegT ======================================
do case_num = 0, case_sum -1
  print("case_num = " +case_num)
  do year = year_s,year_e
    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc" ; outputs/
    f1 = addfile (filename1,"r")

    VegT_tmp   = f1->VegT(0:60*8,{A:B},{C:D})
    do event_num = 0,event_sum-1
      do day = HW_day_s(event_num),HW_day_e(event_num)
        ts = (day-1)*8+4 ;4 - 12pm
        VegT(case_num, 0, day-1, :, :) = VegT_tmp(ts,:,:)
        ts = (day-1)*8+5 ;5 - 12pm
        VegT(case_num, 1, day-1, :, :) = VegT_tmp(ts,:,:) ;dim_avg_n_Wrap(VegT_tmp(ts:ts+1,:,:),0)
      end do
    end do
  end do ; year
end do ; case_num
delete(PFT_mask)

VegT = VegT - 273.15

VegT!0      = "case"
VegT!1      = "time_step"
VegT!2      = "time"
VegT!3      = "lat"
VegT!4      = "lon"
VegT&lat    = lat
VegT&lon    = lon

; ================================== LST =======================================
path_LST = "/g/data/w35/mm3972/data/MODIS_LST/";
Terra    = "MOD11A1.006_1km_aid0001.nc"
Aqua     = "MYD11A1.006_1km_aid0001.nc"

filename2 = path_LST+Terra
f2  = addfile (filename2,"r")

filename3 = path_LST+Aqua
f3  = addfile (filename3,"r")

; Terra passes over Aus between 12 pm and 2 pm

LST_Terra = tofloat(f2->LST_Day_1km(:,{A:B},{C:D}))
LST_Terra = where(LST_Terra .ne. 0, LST_Terra*0.02-273.15,-9999.)
LST_Terra@_FillValue = -9999.

LST_Terra_regrid = linint2_Wrap(doubletofloat(f2->lon({C:D})),doubletofloat(f2->lat({A:B})),LST_Terra,False,lon,lat,0)
LST_Terra_regrid!0 = "time"
LST_Terra_regrid!1 = "lat"
LST_Terra_regrid!2 = "lon"
LST_Terra_regrid&lat = lat
LST_Terra_regrid&lon = lon
LST_Terra_regrid@_FillValue = -9999.
printVarSummary(LST_Terra_regrid)


; Aqua passes over Aus between 3:30 pm and 5:30 pm
LST_Aqua = tofloat(f3->LST_Day_1km(:,{A:B},{C:D}))
LST_Aqua = where(LST_Aqua .ne. 0, LST_Aqua*0.02-273.15,-9999.)
LST_Aqua@_FillValue = -9999.

LST_Aqua_regrid = linint2_Wrap(doubletofloat(f3->lon({C:D})),doubletofloat(f3->lat({A:B})),LST_Aqua,False,lon,lat,0)
LST_Aqua_regrid!0 = "time"
LST_Aqua_regrid!1 = "lat"
LST_Aqua_regrid!2 = "lon"
LST_Aqua_regrid&lat = lat
LST_Aqua_regrid&lon = lon
LST_Aqua_regrid@_FillValue = -9999.
printVarSummary(LST_Aqua_regrid)

; ================================== Tair ======================================
if diff .eq. True then

  Tair     = new((/2,31,dimsizes(lat),dimsizes(lon)/),float,default_val)
  Tair     = 0.

  filename1 = path(0)+"/cable_out_2019_SE_Aus.nc"
  f1 = addfile (filename1,"r")
  Tair_tmp   = f1->Tair(0:60*8,{A:B},{C:D})
  do event_num = 0,event_sum-1
    do day = HW_day_s(event_num),HW_day_e(event_num)
      ts = (day-1)*8+4 ;4 - 12pm
      Tair(0,day-1,:,:) =  Tair_tmp(ts,:,:)
      ts = (day-1)*8+5 ;5 - 3pm
      Tair(1,day-1,:,:) =  Tair_tmp(ts,:,:);dim_avg_n_Wrap(Tair_tmp(ts:ts+1,:,:),0)
    end do
  end do
  delete(Tair_tmp)
  Tair = Tair - 273.15

  VegT(0,:,:,:,:) = VegT(0,:,:,:,:) - Tair
  VegT(1,:,:,:,:) = VegT(1,:,:,:,:) - Tair
  VegT(2,:,:,:,:) = VegT(2,:,:,:,:) - Tair
  LST_Terra_regrid= LST_Terra_regrid - Tair(0,:,:,:)
  LST_Aqua_regrid = LST_Aqua_regrid  - Tair(1,:,:,:)
end if

; ================================= Plotting ===================================
pic = "./plots/"+pic_name
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


plots = new(8,graphic)

if diff .eq. True then
  res@cnLevels      = fspan(-20,20, 41) ;deltaT_CTL_lvl
  res@cnFillPalette = "testcmap";"ViBlGrWhYeOrRe";"WhiteBlueGreenYellowRed"
else
  res@cnLevels      = fspan( 20, 60, 41) ;deltaT_CTL_lvl
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
end if
do day = 0,30
  print("VegT")
  plots(0) = gsn_csm_contour_map(wks,VegT(0,0,day,:,:),res)
  plots(1) = gsn_csm_contour_map(wks,VegT(1,0,day,:,:),res)
  plots(2) = gsn_csm_contour_map(wks,VegT(2,0,day,:,:),res)
  print("LST")
  plots(3) = gsn_csm_contour_map(wks,LST_Terra_regrid(day,:,:),res)

  print("VegT")
  plots(4) = gsn_csm_contour_map(wks,VegT(0,1,day,:,:),res)
  plots(5) = gsn_csm_contour_map(wks,VegT(1,1,day,:,:),res)
  plots(6) = gsn_csm_contour_map(wks,VegT(2,1,day,:,:),res)
  print("LST")
  plots(7) = gsn_csm_contour_map(wks,LST_Aqua_regrid(day,:,:),res)

  do cnt = 0, 7
    txid1   = gsn_create_text(wks, order_name(cnt), txres)
    annoid1 = gsn_add_annotation(plots(cnt), txid1, amres)
  end do
  gsn_panel(wks,(/plots/),(/2,4/),pres)

end do

;
; gsn_panel(wks,(/plots/),(/4,4/),pres)

end


; ; 3: TVeg
; res@cnLevels      = fspan(0,2,21)
; res@cnFillPalette = "WhiteBlueGreenYellowRed";"MPL_RdBu";"NOC_ndvi";"MPL_RdYlGn";"WhiteBlueGreenYellowRed"
; plots(0+count) = gsn_csm_contour_map(wks,TVeg(0,:,:),res)
; delete(res@cnLevels)
;
; res@cnFillPalette = "BlAqGrWh2YeOrReVi22";"ViBlGrWhYeOrRe"
; res@cnLevels      = fspan(-1.,1.,21)
; tmp  = TVeg(1,:,:)
; tmp  = tmp - TVeg(0,:,:)
; plots(1+count) = gsn_csm_contour_map(wks,tmp,res)
; delete(tmp)
;
; tmp  = TVeg(2,:,:)
; tmp  = tmp - TVeg(0,:,:)
; plots(2+count) = gsn_csm_contour_map(wks,tmp,res)
; delete(tmp)
;
; tmp  = TVeg(3,:,:)
; tmp  = tmp - TVeg(0,:,:)
; plots(3+count) = gsn_csm_contour_map(wks,tmp,res)
; delete(tmp)
; delete(res@cnLevels)
;
; count = 12
